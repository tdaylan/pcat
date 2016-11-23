# common imports
from __init__ import *

def retr_psfnwdth(gdat, psfn, frac):
    
    if len(psfn.shape) == 4:
        varioaxi = True
        numboaxi = psfn.shape[3]
        wdth = zeros((gdat.numbener, gdat.numbevtt, numboaxi))
    else:
        varioaxi = False
        numboaxi = psfn.shape[2]
        wdth = zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            for p in arange(numboaxi):
                if varioaxi:
                    psfntemp = psfn[i, :, m, p]
                else:
                    if p > 0:
                        break
                    psfntemp = psfn[i, :, m]
                indxanglgood = argsort(psfntemp)
                intpwdth = max(frac * amax(psfntemp), amin(psfntemp))
                if intpwdth > amin(psfntemp[indxanglgood]) and intpwdth < amax(psfntemp[indxanglgood]):
                    wdthtemp = interp1d(psfntemp[indxanglgood], gdat.binsangl[indxanglgood])(intpwdth)
                if varioaxi:
                    wdth[i, m, p] = wdthtemp
                else:
                    wdth[i, m] = wdthtemp
                        
    return wdth


def retr_spec(gdat, flux, spep=None, spectype=None):

    if isscalar(flux):
        flux = array([flux])

    if spep.ndim == 1:
        spep = spep[None, :]

    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if spectype == 'powr':
            spec = flux[None, :] * gdat.enernorm[:, None]**(-spep[None, :, 0])
        if spectype == 'curv':
            spec = flux[None, :] * gdat.enernorm[:, None]**(-spep[None, :, 0] - gdat.factlogtenerpivt[:, None] * spep[None, :, 1])
        if spectype == 'expo':
            spec = flux[None, :] * gdat.enernorm[:, None]**(-spep[None, :, 0]) * exp(gdat.enerexpofact[:, None] / spep[None, :, 1])

    return spec

# temp -- make a faster version only for indxsampcompcolr
def retr_indx(gdat, indxpntsfull):    

    indxsamplgal = []
    indxsampbgal = []
    indxsampspec = []
    indxsampspep = []
    indxsampcompcolr = []
    for l in gdat.indxpopl:
        indxsamplgaltemp = gdat.indxsampcomp[0] + gdat.maxmnumbcompcuml[l] + array(indxpntsfull[l], dtype=int) * gdat.numbcomp[l]
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], gdat.numbener, 0) + repeat(gdat.indxener, len(indxpntsfull[l])).reshape(gdat.numbener, -1))
        if gdat.numbener > 1:
            indxsampspep.append(indxsamplgaltemp[:, None] + 2 + gdat.numbener + gdat.indxspep[l][None, :])
        indxsampcompcolr.append(repeat(indxsamplgaltemp, gdat.numbcompcolr[l]) + tile(gdat.indxcompcolr[l], len(indxpntsfull[l])))
    
    return indxsamplgal, indxsampbgal, indxsampspec, indxsampspep, indxsampcompcolr


def retr_pntsflux(gdat, lgal, bgal, spec, psfp, psfntype, varioaxi, evalcirc):
   
    numbpnts = lgal.size
    pntsfluxsing = zeros((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for k in range(numbpnts):
        
        if evalcirc:
            indxfluxproxtemp = digitize(spec[gdat.indxenerfluxdist[0], k], gdat.binsfluxprox) - 1
            indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
            indxpixltemp = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
        else:
            indxpixltemp = gdat.indxpixl
        
        # calculate the distance to all pixels from each point source
        dist = retr_angldistunit(gdat, lgal[k], bgal[k], indxpixltemp)

        # evaluate the PSF
        psfn = retr_psfn(gdat, psfp, gdat.indxener, dist, psfntype, gdat.binsoaxi, varioaxi)
            
        # temp
        #psfnintp = interp1d(gdat.binsangl, psfn, axis=1)

        for i in gdat.indxener:
            for m in gdat.indxevtt:
                
                if varioaxi:
                    indxoaxitemp = retr_indxoaxipnts(gdat, lgal[k], bgal[k])
                    psfntemp = psfn[i, :, m, indxoaxitemp]
                else:
                    psfntemp = psfn[i, :, m]
                
                pntsfluxsing[k, i, indxpixltemp, m] = spec[i, k] * psfntemp
                
    # sum contributions from all PS
    pntsflux = sum(pntsfluxsing, 0) 
    
    return pntsflux


def retr_rofi_flux(gdat, bacp, pntsflux, tempindx):
    
    modlflux = pntsflux[tempindx]
    for c in gdat.indxback:
        modlflux += bacp[c, :, None, None] * gdat.backflux[c][tempindx]        
    
    return modlflux


def cdfn_flux_brok(flux, minflux, maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):

    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    fluxunit = norm / (1. - fluxdistsloplowr) * fluxdistbrek**fluxdistsloplowr * (flux**(1. - fluxdistsloplowr) - minmflux**(1. - fluxdistsloplowr))
    indxflux = where(flux >= fluxdistbrek)[0]
    
    if indxflux.size > 0:
        temp = norm * fluxdistbrek**fluxdistsloplowr / (1. - fluxdistsloplowr) * (fluxdistbrek**(1. - fluxdistsloplowr) - minmflux**(1. - fluxdistsloplowr))
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


def icdf_flux_brok(fluxunit, minmflux, maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):
   
    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    fluxunitbrek = norm / (1. - fluxdistsloplowr) * fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - minmflux**(1. - fluxdistsloplowr))
    flux = (fluxunit * (1. - fluxdistsloplowr) / norm / fluxdistbrek**fluxdistsloplowr + minmflux**(1. - fluxdistsloplowr))**(1. / (1. - fluxdistsloplowr))
    indxfluxunit = where(fluxunit >= fluxunitbrek)[0]
    
    if indxfluxunit.size > 0:
        temp = norm * fluxdistbrek**fluxdistsloplowr / (1. - fluxdistsloplowr) * (fluxdistbrek**(1. - fluxdistsloplowr) - minmflux**(1. - fluxdistsloplowr))
        flux[indxfluxunit] = ((fluxunit[indxfluxunit] - temp) * (1. - fluxdistslopuppr) / norm / fluxdistbrek**fluxdistslopuppr + \
                                                                                            fluxdistbrek**(1. - fluxdistslopuppr))**(1. / (1. - fluxdistslopuppr))

    return flux


def cdfn_flux_powr(flux, minmflux, maxmflux, fluxdistslop):
        
    fluxunit = (flux**(1. - fluxdistslop) - minmflux**(1. - fluxdistslop)) / (maxmflux**(1. - fluxdistslop) - minmflux**(1. - fluxdistslop))
        
    return fluxunit


def icdf_flux_powr(fluxunit, minmflux, maxmflux, fluxdistslop):

    flux = (fluxunit * (maxmflux**(1. - fluxdistslop) - minmflux**(1. - fluxdistslop)) + minmflux**(1. - fluxdistslop))**(1. / (1. - fluxdistslop))

    return flux


def pdfn_flux_powr(gdat, flux, fluxdistslop):
  
    norm = (1. - fluxdistslop) / (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop))
    
    pdfn = norm * flux**(-fluxdistslop)
    
    return pdfn


def cdfn_self(para, minmpara, factpara):
    
    paraunit = (para - minmpara) / factpara
    
    return paraunit


def icdf_self(paraunit, minmpara, factpara):
    
    para = factpara * paraunit + minmpara
    
    return para


def cdfn_gaus(para, meanpara, stdvpara):
   
    paraunit = 0.5  * (1. + sp.special.erf((para - meanpara) / sqrt(2) / stdvpara))
    
    return paraunit


def icdf_gaus(paraunit, meanpara, stdvpara):
    
    para = meanpara + stdvpara * sqrt(2) * sp.special.erfinv(2. * paraunit - 1.)

    return para


def cdfn_eerr(para, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    tranpara = (para - meanpara) / stdvpara
    cdfnnormpara = 0.5 * (sp.special.erf(tranpara / sqrt(2.)) + 1.)
    paraunit = (cdfnnormpara - cdfnnormminm) / cdfnnormdiff

    return paraunit


def icdf_eerr(paraunit, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    cdfnnormpara = paraunit * cdfnnormdiff + cdfnnormminm
    tranpara = sp.special.erfinv(2. * cdfnnormpara - 1.) * sqrt(2)
    para = tranpara * stdvpara + meanpara
   
    return para


def cdfn_logt(para, minmpara, factpara):

    paraunit = log(para / minmpara) / factpara

    return paraunit


def icdf_logt(paraunit, minmpara, factpara):

    para = exp(paraunit * factpara) * minmpara

    return para


def cdfn_atan(para, minmpara, factpara):
    
    paraunit = (arctan(para) - arctan(minmpara)) / factpara
    
    return paraunit


def icdf_atan(paraunit, minmpara, factpara):

    para = tan(factpara * paraunit + arctan(minmpara))
    
    return para


def cdfn_fixp(gdat, fixp, thisindxfixp):

    scalfixp = gdat.scalfixp[thisindxfixp]
    if scalfixp == 'self' or scalfixp == 'logt' or scalfixp == 'atan':
        minmfixp = gdat.minmfixp[thisindxfixp]
        factfixp = gdat.factfixp[thisindxfixp]
        if scalfixp == 'self':
            fixpunit = cdfn_self(fixp, minmfixp, factfixp)
        elif scalfixp == 'logt':
            fixpunit = cdfn_logt(fixp, minmfixp, factfixp)
        elif scalfixp == 'atan':
            fixpunit = cdfn_atan(fixp, minmfixp, factfixp)
    elif scalfixp == 'gaus' or scalfixp == 'eerr':
        meanfixp = gdat.meanfixp[thisindxfixp]
        stdvfixp = gdat.stdvfixp[thisindxfixp]
        if scalfixp == 'eerr':
            cdfnminmfixp = gdat.cdfnminmfixp[thisindxfixp]
            cdfndifffixp = gdat.cdfndifffixp[thisindxfixp]
            fixpunit = cdfn_eerr(fixp, meanfixp, stdvfixp, cdfnminmfixp, cdfndifffixp)
        else:
            fixpunit = cdfn_gaus(fixp, meanfixp, stdvfixp)
    elif scalfixp == 'pois':
        fixpunit = fixp

    return fixpunit


def icdf_fixp(gdat, fixpunit, thisindxfixp):

    scalfixp = gdat.scalfixp[thisindxfixp]
    if scalfixp == 'self' or scalfixp == 'logt' or scalfixp == 'atan':
        minmfixp = gdat.minmfixp[thisindxfixp]
        factfixp = gdat.factfixp[thisindxfixp]
        if scalfixp == 'self':
            fixp = icdf_self(fixpunit, minmfixp, factfixp)
        elif scalfixp == 'logt':
            fixp = icdf_logt(fixpunit, minmfixp, factfixp)
        elif scalfixp == 'atan':
            fixp = icdf_atan(fixpunit, minmfixp, factfixp)
    elif scalfixp == 'gaus' or scalfixp == 'eerr':
        meanfixp = gdat.meanfixp[thisindxfixp]
        stdvfixp = gdat.stdvfixp[thisindxfixp]
        if scalfixp == 'eerr':
            cdfnminmfixp = gdat.cdfnminmfixp[thisindxfixp]
            cdfndifffixp = gdat.cdfndifffixp[thisindxfixp]
            fixp = icdf_eerr(fixpunit, meanfixp, stdvfixp, cdfnminmfixp, cdfndifffixp)
        else:
            fixp = icdf_gaus(fixpunit, meanfixp, stdvfixp)
    elif scalfixp == 'pois':
        fixp = fixpunit
    else:
        raise Exception('Scaling of the parameter is unrecognized.')

    return fixp


def retr_thisindxprop(gdat, gdatmodi):

    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal,  gdatmodi.thisindxsampspec, gdatmodi.thisindxsampspep, \
                                                                            gdatmodi.thisindxsampcompcolr = retr_indx(gdat, gdatmodi.thisindxpntsfull)
 
    gdatmodi.propcomp = False
    if rand() < gdat.probproptran:
        gdatmodi.indxpoplmodi = choice(gdat.indxpopl)
        if rand() < gdat.probpropbrth:
            if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropdeth
            elif gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.minmnumbpnts[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropbrth
            else:
                gdatmodi.thisindxprop = choice(gdat.indxpropbrde, p=gdat.probtran)
        else:
            if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropmerg
            elif gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.minmnumbpnts[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropbrth
            else:
                gdatmodi.thisindxprop = choice(gdat.indxpropspme, p=gdat.probtran)
            
        gdatmodi.proptran = True
        gdatmodi.propfixp = False
    else:
        indxsampfull = concatenate((concatenate(gdatmodi.thisindxsampcompcolr), gdat.indxfixpusedprim))
        gdatmodi.indxsampmodi = choice(indxsampfull)
        gdatmodi.proptran = False
        
        if gdatmodi.indxsampmodi < gdat.indxsampcomp[0]:
            gdatmodi.propfixp = True
        else:
            gdatmodi.propfixp = False
    
        if gdatmodi.indxsampmodi in gdat.indxfixpmeanpnts:
            gdatmodi.thisindxprop = gdat.indxpropmeanpnts
        elif gdatmodi.indxsampmodi in gdat.indxfixpfluxdist:
            gdatmodi.thisindxprop = gdat.indxpropfluxdist
        elif gdatmodi.indxsampmodi in gdat.indxfixppsfp:
            gdatmodi.thisindxprop = gdat.indxproppsfp
        elif gdatmodi.indxsampmodi in gdat.indxfixpbacp:
            gdatmodi.thisindxprop = gdat.indxpropbacp
        elif gdatmodi.indxsampmodi in gdat.indxfixplenp:
            gdatmodi.thisindxprop = gdat.indxproplenp
        else:
            indxsampdiff = gdatmodi.indxsampmodi - gdat.indxsampcomp[0]
            gdatmodi.indxpoplmodi = argmin(where(indxsampdiff // gdat.maxmnumbcompcumr == 0)[0])
            indxsamptemp = indxsampdiff - gdat.maxmnumbcompcuml[gdatmodi.indxpoplmodi]
            gdatmodi.indxpntsmodi = indxsamptemp // gdat.numbcomp[gdatmodi.indxpoplmodi]
            gdatmodi.indxpntsfullmodi = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].index(gdatmodi.indxpntsmodi)
            gdatmodi.indxcompmodi = indxsamptemp % gdat.numbcomp[gdatmodi.indxpoplmodi]

            # temp
            if gdatmodi.indxsampmodi in gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxproplgal
            elif gdatmodi.indxsampmodi in gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropbgal
            elif gdatmodi.indxsampmodi in gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropflux
            elif gdatmodi.indxsampmodi in gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi]:
                gdatmodi.thisindxprop = gdat.indxpropspep
                gdatmodi.indxspepmodi = gdatmodi.indxcompmodi % gdat.indxcompspep[0]
            gdatmodi.propcomp = True
    
    if gdat.verbtype > 1:
        print 
        print 'retr_thisindxprop()'
        print 'gdatmodi.proptran'
        print gdatmodi.proptran
        print 'propfixp'
        print gdatmodi.propfixp
        print 'gdatmodi.thisindxprop'
        print gdatmodi.thisindxprop
        print 'gdat.strgprop[gdatmodi.thisindxprop]'
        print gdat.strgprop[gdatmodi.thisindxprop]
        if gdatmodi.proptran and gdat.verbtype > 1:
            print 'gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]'
            print gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]
            print 'indxpoplmodi'
            print gdatmodi.indxpoplmodi
        
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
        if gdat.pixltype == 'unbd':
            gdatmodi.thisllik = gdat.numbdatasamp * log(gdatmodi.thismodlfluxtotl) - gdatmodi.thismodlfluxtotl + log(gdatmodi.thismodlflux)
        else:
            if gdat.liketype == 'pois':
                gdatmodi.thisllik = gdat.datacnts * log(gdatmodi.thismodlcnts) - gdatmodi.thismodlcnts
            if gdat.liketype == 'gaus':
                gdatmodi.thisllik = -0.5 * (gdat.datacnts - gdatmodi.thismodlcnts)**2 / gdat.datacnts
  
    elif gdatmodi.thisindxprop in gdat.indxpropllik:

        timeinit = gdat.functime()
        
        # PSF parameter change 
        if gdatmodi.thisindxprop == gdat.indxproppsfp and gdat.evalpsfnpnts:
            
            gdatmodi.numbpntsmodi = int(sum(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts]))
            # construct the proposed PSF
            gdatmodi.nextpsfn = retr_psfn(gdat, gdatmodi.nextsampvarb[gdat.indxfixppsfp], gdat.indxener, gdat.binsangl, gdat.modlpsfntype, gdat.binsoaxi, gdat.modlvarioaxi)
            if gdat.modlvarioaxi:
                for k in gdat.indxoaxi:
                    gdatmodi.nextpsfnintp[k] = interp1d(gdat.binsangl, gdatmodi.nextpsfn[:, :, :, k], axis=1)
            else:
                gdatmodi.nextpsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
            

        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 0] = timefinl - timeinit
            
        if gdat.verbtype > 1:
            print 'gdatmodi.numbpntsmodi'
            print gdatmodi.numbpntsmodi
            print

        # load convenience variables
        timeinit = gdat.functime()
        if gdatmodi.thisindxprop == gdat.indxproppsfp:
            lgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsamplgal)]
            bgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampbgal)]
            spec = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdatmodi.indxenermodi, :]]
        if gdatmodi.thisindxprop >= gdat.indxpropbrth:
            lgal = gdatmodi.modilgal[:gdatmodi.numbpntsmodi]
            bgal = gdatmodi.modibgal[:gdatmodi.numbpntsmodi]
            spec = gdatmodi.modispec[meshgrid(gdatmodi.indxenermodi, arange(gdatmodi.numbpntsmodi), indexing='ij')]
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 0] = timefinl - timeinit
        
        # determine pixels over which to evaluate the log-likelihood
        timeinit = gdat.functime()
        if gdatmodi.thisindxprop == gdat.indxproplenp or gdatmodi.thisindxprop == gdat.indxpropbacp or not gdat.evalcirc:
            gdatmodi.indxpixlmodi = gdat.indxpixl
        else:
            thisindxpixlprox = []
            for k in range(gdatmodi.numbpntsmodi):
                # temp -- this may not work for extreme color PS!
                # take the flux at the pivot energy
                if gdatmodi.thisindxprop == gdat.indxproppsfp:
                    fluxtemp = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdat.indxenerfluxdist, k]]
                else:
                    fluxtemp = gdatmodi.modispec[gdat.indxenerfluxdist, k]
                
                # find the flux index
                # temp
                indxfluxproxtemp = digitize(fabs(fluxtemp), gdat.binsfluxprox) - 1
                indxpixltemp = retr_indxpixl(gdat, bgal[k], lgal[k])
                thisindxpixlprox.append(gdat.indxpixlprox[indxfluxproxtemp][indxpixltemp])
            gdatmodi.indxpixlmodi = unique(concatenate(thisindxpixlprox))
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 1] = timefinl - timeinit
        
        # construct the mesh grid for likelihood evaluation
        timeinit = gdat.functime()
        
        if not gdat.evalcirc and (gdat.numbener == 1 or gdatmodi.propcomp):
            gdatmodi.indxcubemodi = gdat.indxcube
        else:
            gdatmodi.indxcubemodi = meshgrid(gdatmodi.indxenermodi, gdatmodi.indxpixlmodi, gdat.indxevtt, indexing='ij')
        timefinl = gdat.functime()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 2] = timefinl - timeinit

        if gdat.pntstype == 'lght':
            # update the model PS flux map, if needed
            timeinit = gdat.functime()
            if gdatmodi.thisindxprop == gdat.indxproppsfp or gdatmodi.thisindxprop >= gdat.indxpropbrth:
                
                # copy the previous PS flux map
                gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] = copy(gdatmodi.thispntsflux[gdatmodi.indxcubemodi])
                    
                ## when evaluating the PSF, avoid copying the current PS fluxes unnecessarily
                if gdatmodi.thisindxprop == gdat.indxproppsfp:
                    numbrept = 2
                else:
                    numbrept = 1

                for n in range(numbrept):
                    
                    # grab the PSF interpolating function
                    if gdatmodi.thisindxprop == gdat.indxproppsfp:
                        if n == 0:
                            psfnintp = gdatmodi.thispsfnintp
                        else:
                            psfnintp = gdatmodi.nextpsfnintp
                    else:
                        psfnintp = gdatmodi.thispsfnintp

                    for k in range(gdatmodi.numbpntsmodi):
                        
                        # calculate the distance to the pixels to be updated
                        if gdat.evalcirc:
                            dist = retr_angldistunit(gdat, lgal[k], bgal[k], thisindxpixlprox[k])
                        else:
                            dist = retr_angldistunit(gdat, lgal[k], bgal[k], gdat.indxpixl)
                            
                        if gdat.verbtype > 1:
                            print 'dist'
                            print amin(dist) * gdat.anglfact, amax(dist) * gdat.anglfact
                            print 
                        
                        # interpolate the PSF for each PS over the set of data pixels to be updated
                        if gdat.modlvarioaxi:
                            indxoaxitemp = retr_indxoaxipnts(gdat, lgal[k], bgal[k])
                            psfn = psfnintp[indxoaxitemp](dist)
                        else:
                            psfn = psfnintp(dist)
                        
                        for i in range(gdatmodi.indxenermodi.size):
    
                            # expedite PSF proposals
                            if gdatmodi.thisindxprop == gdat.indxproppsfp:
                                if n == 0:
                                    spectemp = -spec[i, k]
                                else:
                                    spectemp = spec[i, k]
                            else:
                                spectemp = spec[i, k]

                            # add the contribution of the PS to the the proposed flux map
                            if gdat.evalcirc:
                                gdatmodi.nextpntsflux[gdatmodi.indxenermodi[i], thisindxpixlprox[k], :] += spectemp * psfn[gdatmodi.indxenermodi[i], :, :]
                            else:
                                gdatmodi.nextpntsflux[gdatmodi.indxenermodi[i], :, :] += spectemp * psfn[gdatmodi.indxenermodi[i], :, :]
                                
            timefinl = gdat.functime()
            gdatmodi.listchrollik[gdatmodi.cntrswep, 3] = timefinl - timeinit
            
            # update the total model flux map
            timeinit = gdat.functime()
            indxtemp = meshgrid(gdat.indxback, gdatmodi.indxenermodi, indexing='ij')
            if gdatmodi.thisindxprop == gdat.indxpropbacp:
                bacp = gdatmodi.nextsampvarb[gdat.indxfixpbacp[indxtemp]]
                pntsflux = gdatmodi.thispntsflux
            if gdatmodi.thisindxprop == gdat.indxproppsfp or gdatmodi.thisindxprop >= gdat.indxpropbrth:
                bacp = gdatmodi.thissampvarb[gdat.indxfixpbacp[indxtemp]]
                pntsflux = gdatmodi.nextpntsflux
            gdatmodi.nextmodlflux[gdatmodi.indxcubemodi] = retr_rofi_flux(gdat, bacp, pntsflux, gdatmodi.indxcubemodi)
            
            timefinl = gdat.functime()
            gdatmodi.listchrollik[gdatmodi.cntrswep, 4] = timefinl - timeinit
    
        if gdatmodi.thisindxprop in gdat.indxfixplenp:
            
            timeinit = gdat.functime()
            
            # calculate the deflection matrix
            gdatmodi.nextdefl = copy(gdatmodi.thisdefl)
            gdatmodi.indxcubemodi = gdat.indxcube
            for k in range(gdatmodi.numbpntsmodi):
                if gdatmodi.modispec[0, k] < 0:
                    modispectemp = abs(gdatmodi.modispec[0, k])
                    facttemp = -1.
                else:
                    modispectemp = gdatmodi.modispec[0, k]
                    facttemp = 1.
                
                gdatmodi.thislensmodl = franlens.LensModel(gdat.mocklenstype, gdatmodi.modilgal[k], gdatmodi.modibgal[k], gdat.mockellplens, \
                                                                                         gdat.mockangllens, gdat.mocksherlens, gdat.mocksanglens, modispectemp)
                
                gdatmodi.nextdefl += facttemp * gdatmodi.thislensmodl.deflection(gdat.lgalgrid, gdat.bgalgrid)
           
            # calculate the lensed image
            if gdatmodi.thisindxprop in gdat.indxfixpsour:
                sourobjt = franlens.Source(gdat.mocksourtype, gdatmodi.nextsampvarb[gdat.indxfixplgalsour], \
                                                              gdatmodi.nextsampvarb[gdat.indxfixpbgalsour], \
                                                              gdatmodi.nextsampvarb[gdat.indxfixpfluxsour], \
                                                              gdatmodi.nextsampvarb[gdat.indxfixpsizesour], \
                                                              gdatmodi.nextsampvarb[gdat.indxfixpratisour], \
                                                              gdatmodi.nextsampvarb[gdat.indxfixpanglsour])
            else:
                sourobjt = gdatmodi.thissourobjt
            gdatmodi.nextpntsflux[0, :, 0] = sourobjt.brightness(gdat.lgalgrid - gdatmodi.nextdefl[:, 0], gdat.bgalgrid - gdatmodi.nextdefl[:, 1])

        # convolve the lensed image with the PSF
        if gdat.pntstype == 'lens':
            if gdatmodi.thisindxprop in gdat.indxfixppsfp:
                sigc = gdatmodi.thissampvarb[gdat.indxfixppsfp[0]]
            else:
                sigc = gdatmodi.nextsampvarb[gdat.indxfixppsfp[0]]
            gdatmodi.nextpntsflux[0, :, 0] = convolve(gdatmodi.nextpntsflux[0, :, 0], sigc)
                
            gdatmodi.nextmodlflux = gdatmodi.nextpntsflux

            timefinl = gdat.functime()
            gdatmodi.listchrollik[gdatmodi.cntrswep, 0] = timefinl - timeinit
    
        if gdat.correxpo:
            # calculate the count map
            timeinit = gdat.functime()
            gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi] = gdatmodi.nextmodlflux[gdatmodi.indxcubemodi] * gdat.expo[gdatmodi.indxcubemodi] * \
                                                                                                      gdat.apix * gdat.diffener[gdatmodi.indxenermodi, None, None] # [1]
            timefinl = gdat.functime()
            if gdat.pntstype == 'lght':
                gdatmodi.listchrollik[gdatmodi.cntrswep, 5] = timefinl - timeinit
            if gdat.pntstype == 'lens':
                gdatmodi.listchrollik[gdatmodi.cntrswep, 1] = timefinl - timeinit
    
        # temp
        if gdat.pixltype == 'unbd':
            gdatmodi.nextmodlfluxtotl = gdatmodi.thismodlfluxtotl + sum(gdatmodi.modispec[gdatmodi.indxenermodi, :gdatmodi.numbpntsmodi])
        
        # calculate the log-likelihood difference over the modified data cubes
        timeinit = gdat.functime()
        if gdat.pixltype == 'unbd':
            gdatmodi.nextllik = gdat.numbdatasamp * log(gdatmodi.nextmodlfluxtotl) - gdatmodi.nextmodlfluxtotl + log(gdatmodi.nextmodlflux)
        else:
            if gdat.liketype == 'pois':
                gdatmodi.nextllik[gdatmodi.indxcubemodi] = gdat.datacnts[gdatmodi.indxcubemodi] * log(gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi]) \
                                                                                                                                 - gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi]
            if gdat.liketype == 'gaus':
                gdatmodi.nextllik[gdatmodi.indxcubemodi] = -0.5 * (gdat.datacnts[gdatmodi.indxcubemodi] - gdatmodi.nextmodlcnts[gdatmodi.indxcubemodi])**2 / \
                                                                                                                                    gdat.datacnts[gdatmodi.indxcubemodi]
            
            gdatmodi.deltllik = sum(gdatmodi.nextllik[gdatmodi.indxcubemodi] - gdatmodi.thisllik[gdatmodi.indxcubemodi])
        timefinl = gdat.functime()
        if gdat.pntstype == 'lght':
            gdatmodi.listchrollik[gdatmodi.cntrswep, 6] = timefinl - timeinit
        if gdat.pntstype == 'lens':
            gdatmodi.listchrollik[gdatmodi.cntrswep, 2] = timefinl - timeinit
    
        if gdat.diagmode:
            if not isfinite(gdatmodi.nextllik[gdatmodi.indxcubemodi]).any():
                warnings.warn('Log-likelihood went NAN!')

        # temp
        if False:
            temppntsflux, temppntscnts, tempmodlflux, tempmodlcnts = retr_maps(gdat, list(gdatmodi.thisindxpntsfull), copy(gdatmodi.thissampvarb))
            gdatmodi.thiserrrcnts = gdatmodi.thispntscnts - temppntscnts
            gdatmodi.thiserrr = zeros_like(gdatmodi.thiserrrcnts)
            indxcubegood = where(temppntscnts > 0.)
            gdatmodi.thiserrr[indxcubegood] = 100. * gdatmodi.thiserrrcnts[indxcubegood] / temppntscnts[indxcubegood]
            
            facttemp = gdat.expo[gdatmodi.indxcubemodi] * gdat.diffener[gdatmodi.indxenermodi] * gdat.apix
    
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = temppntscnts[gdatmodi.indxcubemodi]
            path = gdat.pathdiag + '0temppntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, pixltype=gdat.pixltype, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, titl=gdat.strgprop[gdatmodi.thisindxprop], \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
            
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = gdatmodi.thiserrrcnts[gdatmodi.indxcubemodi]
            path = gdat.pathdiag + '1errrcnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, titl=gdat.strgprop[gdatmodi.thisindxprop], \
                                                                              resi=True, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
    
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = gdatmodi.thiserrr[gdatmodi.indxcubemodi]
            path = gdat.pathdiag + '2errr_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, titl=gdat.strgprop[gdatmodi.thisindxprop], \
                                                                              resi=True, \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
    
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = gdatmodi.thispntsflux[gdatmodi.indxcubemodi] * facttemp
            path = gdat.pathdiag + '3thispntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, titl=gdat.strgprop[gdatmodi.thisindxprop], \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
        
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] * facttemp
            path = gdat.pathdiag + '4nextpntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, titl=gdat.strgprop[gdatmodi.thisindxprop], \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
        
            temp = zeros(gdat.numbpixl)
            temp[gdatmodi.indxpixlmodi] = (gdatmodi.nextpntsflux[gdatmodi.indxcubemodi] - gdatmodi.thispntsflux[gdatmodi.indxcubemodi]) * facttemp
            path = gdat.pathdiag + '5diffpntscnts_%09d.pdf' % gdatmodi.cntrswep
            tdpy.util.plot_maps(path, temp, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, titl=gdat.strgprop[gdatmodi.thisindxprop], \
                                                                              minmlgal=0.95*gdat.anglfact*gdat.minmlgal, maxmlgal=0.95*gdat.anglfact*gdat.maxmlgal, \
                                                                              minmbgal=0.95*gdat.anglfact*gdat.minmbgal, maxmbgal=0.95*gdat.anglfact*gdat.maxmbgal)
            
            if amax(fabs(gdatmodi.thiserrrcnts)) > 0.1:
                print 
                raise Exception('Approximation error in calculating the PS flux map is above the tolerance level.')
            if gdatmodi.cntrswep > 50:
                raise Exception('')
                
    else:
        gdatmodi.deltllik = 0.
        

def retr_cntsbackfwhm(gdat, bacp, fwhm):

    varioaxi = len(fwhm.shape) == 3
    cntsbackfwhm = zeros_like(fwhm)
    for c in gdat.indxback:
        if varioaxi:
            cntsback = bacp[c, :, None, None, None] * gdat.backflux[c][:, :, :, None] * gdat.expo[:, :, :, None] * \
                                                                                                gdat.diffener[:, None, None, None] * pi * fwhm[:, None, :, :]**2 / 4.
        else:
            cntsback = bacp[c, :, None, None] * gdat.backflux[c] * gdat.expo * gdat.diffener[:, None, None] * pi * fwhm[:, None, :]**2 / 4.
        cntsbackfwhm += mean(cntsback, 1)

    return cntsbackfwhm


def retr_sigm(gdat, cnts, cntsbackfwhm, lgal=None, bgal=None):
   
    varioaxi = len(cntsbackfwhm.shape) == 3
    if cnts.ndim == 2:
        if varioaxi:
            sigm = cnts / sum(cntsbackfwhm[:, :, 0], 1)[:, None]
        else:
            sigm = cnts / sum(cntsbackfwhm, 1)[:, None]
    else:
        if varioaxi:
            indxoaxitemp = retr_indxoaxipnts(gdat, lgal, bgal)
            sigm = cnts / swapaxes(cntsbackfwhm[:, :, indxoaxitemp], 1, 2)
        else:
            sigm = cnts / cntsbackfwhm[:, None, :]

    return sigm


def retr_fluxlpribind(gdatmodi, gdat, l):
    
    if gdat.fluxdisttype[l] == 'powr':
        fluxhistmodl = gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[l]] * gdat.diffflux * pdfn_flux_powr(gdat, gdat.meanflux, \
                                                                               gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]])
    if gdat.fluxdisttype[l] == 'brok':
        fluxhistmodl = gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[l]] * gdat.diffflux * pdfn_flux_brok(gdat, gdat.meanflux, \
                                                                               gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[l]], \
                                                                               gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[l]], \
                                                                               gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[l]])
    fluxhistdata = histogram(thisflux, gdat.binsflux)[0]
    fluxlpribind = retr_probpois(fluxhistdata, fluxhistmodl)
    
    return fluxlpribind


def retr_fluxlpri(gdatmodi, gdat, l):
    
    if gdat.fluxdisttype[l] == 'powr':
        fluxlpri = sum(log(pdfn_flux_powr(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]], \
                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]])))
    if gdat.fluxdisttype[l] == 'brok':
        fluxlpri = sum(log(pdfn_flux_brok(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]], \
                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[l]], \
                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[l]], \
                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[l]])))
    
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
                gdatmodi.thislpri[l, 0] = retr_probpois(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]], gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[l]])
                gdatmodi.thislpri[l, 1] = retr_fluxlpri(gdatmodi, gdat, l)
        gdatmodi.thislpritotl = sum(gdatmodi.thislpri)
        gdatmodi.nextlpri = copy(gdatmodi.thislpri)
    else:
        
        # initialize the delta log-prior
        gdatmodi.deltlpri = 0.

        # determine which type of delta log-prior is to be calculated
        boolupdtmeanpnts = gdatmodi.thisindxprop == gdat.indxpropmeanpnts
        boolupdtfluxdist = gdatmodi.thisindxprop == gdat.indxpropfluxdist
        boolupdtnumbpnts = gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg

        # calculate contributions to the delta log-prior
        if boolupdtmeanpnts or boolupdtfluxdist or boolupdtnumbpnts:

            # penalty term due to the number of degrees of freedom
            if boolupdtnumbpnts:
                if gdatmodi.thisindxprop == gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxpropsplt:
                    deltdoff = -gdat.numbcompcolr[gdatmodi.indxpoplmodi]
                else:
                    deltdoff = gdat.numbcompcolr[gdatmodi.indxpoplmodi]
                gdatmodi.deltlpri += gdat.priofactdoff * deltdoff

            # binned flux prior
            # temp -- binned prior currently does not work with splits and merges!
            if gdat.bindprio:
               
                # mean number of PS
                if boolupdtmeanpnts:
                    meanpnts = gdatmodi.nextsampvarb[gdat.indxfixpmeanpnts[gdatmodi.indxpoplmodi]]
                else:
                    meanpnts = gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[gdatmodi.indxpoplmodi]]
                
                # hyperparameters on the flux distribution
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                if boolupdtfluxdist:
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                        fluxdistslop = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                            fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]]
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                            fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
                            fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                
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
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] = retr_probpois(gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[gdatmodi.indxpoplmodi]])
                    else:
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] = retr_probpois(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.nextsampvarb[gdat.indxfixpmeanpnts[gdatmodi.indxpoplmodi]])
                    gdatmodi.deltlpri += gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, 0]
                    
                if boolupdtfluxdist:
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] = sum(log(pdfn_flux_powr(gdat, \
                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]], \
                                                                    gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])))
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                        fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]]
                        fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                        if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                            fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]]
                        if gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                            fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        if gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
                            fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] = sum(log(pdfn_flux_brok(gdat, \
                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]], \
                                                    fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)))
                    gdatmodi.deltlpri += gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, 1]
       
            if gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg:
                
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    lprbfrst = log(pdfn_flux_powr(gdat, gdatmodi.fluxfrst, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]))
                    lprbseco = log(pdfn_flux_powr(gdat, gdatmodi.fluxseco, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]))
                    lprbpare = log(pdfn_flux_powr(gdat, gdatmodi.fluxpare, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]))
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    lprbfrst += log(pdfn_flux_brok(gdat, gdatmodi.fluxfrst, \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]))
                    
                    lprbseco += log(pdfn_flux_brok(gdat, gdatmodi.fluxseco, \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]))

                    lprbpare -= log(pdfn_flux_brok(gdat, gdatmodi.fluxpare, \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
                                       gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]))
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
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampspep, indxsampcompcolr = retr_indx(gdat, indxpntsfull)    
    
    sampvarb = zeros_like(samp)
    for k in gdat.indxfixp:
        sampvarb[gdat.indxfixp[k]] = icdf_fixp(gdat, samp[gdat.indxfixp[k]], k)
    
    for l in gdat.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl) 
        if gdat.fluxdisttype[l] == 'powr':
            sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_powr(samp[indxsampspec[l][gdat.indxenerfluxdist, :]], gdat.minmflux, gdat.maxmflux, \
                                                                                                                                            sampvarb[gdat.indxfixpfluxdistslop[l]])
        if gdat.fluxdisttype[l] == 'brok':
            fluxunit = samp[indxsampspec[l][gdat.indxenerfluxdist[0], :]]
            fluxdistbrek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
            fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
            fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
            sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_brok(fluxunit, gdat.minmflux, gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        
        if gdat.numbener > 1:
            sampvarb[indxsampspep[l][:, 0]] = icdf_gaus(samp[indxsampspep[l][:, 0]], gdat.sinddistmean[l], gdat.sinddiststdv[l])
            if gdat.spectype[l] == 'curv':
                sampvarb[indxsampspep[l][:, 1]] = icdf_gaus(samp[indxsampspep[l][:, 1]], gdat.curvddistmean[l], gdat.curvdiststdv[l])
            if gdat.spectype[l] == 'expo':
                sampvarb[indxsampspep[l][:, 1]] = icdf_logt(samp[indxsampspep[l][:, 1]], gdat.minmener, gdat.factener)
        
        if gdat.numbener > 1:   
            sampvarb[indxsampspec[l]] = retr_spec(gdat, sampvarb[indxsampspec[l][gdat.indxenerfluxdist[0], :]], spep=sampvarb[indxsampspep[l]], spectype=gdat.spectype[l])

    return sampvarb
    

def retr_maps(gdat, indxpntsfull, sampvarb, evalcirc):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampspep, indxsampcompcolr = retr_indx(gdat, indxpntsfull)    
     
    listspectemp = []
    for l in gdat.indxpopl:
        listspectemp.append(sampvarb[indxsampspec[l]])

    pntsflux = retr_pntsflux(gdat, sampvarb[concatenate(indxsamplgal)], sampvarb[concatenate(indxsampbgal)], \
                                                            concatenate(listspectemp, axis=1), sampvarb[gdat.indxfixppsfp], gdat.modlpsfntype, gdat.modlvarioaxi, evalcirc)
    
    if gdat.backemis:
        totlflux = retr_rofi_flux(gdat, sampvarb[gdat.indxfixpbacp], pntsflux, gdat.indxcube)
    else:
        totlflux = pntsflux
        
    if gdat.pixltype == 'unbd':
        
        totlcntstotl = totlflux * gdat.apix

        return pntsflux, totlflux, totlcntstotl
    
    else:
    
        pntscnts = pntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
        totlcnts = totlflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
        
        return pntsflux, pntscnts, totlflux, totlcnts


def retr_mrkrsize(gdat, flux):

    mrkrsize = (log(flux) - log(gdat.minmflux)) / (log(gdat.maxmflux) - log(gdat.minmflux)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


def retr_hubbpsfn(gdat):

    gdat.truepsfp = array([0.05, 0.05]) / gdat.anglfact


def retr_sdsspsfn(gdat):
   
    gdat.truepsfp = array([0.25 / gdat.anglfact, 1.7e6, 1.9, 0.25 / gdat.anglfact, 2.1e6, 2.])


def retr_chanpsfn(gdat):

    gdat.truepsfp = array([0.3 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.6e-1, 2.])
   

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

    numbpsfpscal = 3
    numbpsfpform = 5
    
    fermscal = zeros((gdat.numbevtt, numbpsfpscal))
    fermform = zeros((gdat.numbener, gdat.numbevtt, numbpsfpform))
    
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
        for k in range(numbpsfpform):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)
        
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)

    # store the fermi PSF parameters
    gdat.truepsfp = zeros((gdat.numbener * numbpsfpform * gdat.numbevtt))
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            gdat.truepsfp[indxfermpsfptemp] = fermform[:, m, k]

    # calculate the scale factor
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    

def updt_samp(gdat, gdatmodi):
    
    # update the sample vector
    if gdatmodi.thisindxprop != gdat.indxpropdeth:
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
    
    # update the log-prior
    if gdatmodi.thisindxprop in gdat.indxproplpri:
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]
        gdatmodi.thislpritotl = sum(gdatmodi.thislpri)
    
    # update the log-likelihood
    if gdatmodi.thisindxprop in gdat.indxpropllik:
        gdatmodi.thisllik[gdatmodi.indxcubemodi] = gdatmodi.nextllik[gdatmodi.indxcubemodi]
        gdatmodi.thislliktotl = sum(gdatmodi.thisllik)
        
    ## update the unit sample vector for the change in the flux distribution
    if gdatmodi.thisindxprop == gdat.indxpropfluxdist:
        flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]]
        gdatmodi.indxfluxdistpara = (gdatmodi.indxsampmodi - gdat.indxfixpfluxdist[0])
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
            fluxdistslop = gdatmodi.nextsampvarb[gdat.indxfixpfluxdist[gdatmodi.indxfluxdistpara]]
            fluxunit = cdfn_flux_powr(flux, gdat.minmflux, gdat.maxmflux, fluxdistslop)
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
            if gdatmodi.indxfluxdistpara == gdatmodi.indxfluxdistbrek:
                fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistbrek]]
                fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistsloplowr]]
                fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistslopuppr]]
            elif gdatmodi.indxfluxdistpara == gdatmodi.indxfluxdistsloplowr:
                fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistbrek]]
                fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistsloplowr]]
                fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistslopuppr]]
            else:
                fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistbrek]]
                fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistsloplowr]]
                fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxfixpfluxdist[gdatmodi.indxpoplmodi][gdatmodi.indxfluxdistslopuppr]]
            fluxunit = cdfn_flux_brok(flux, gdat.minmflux, gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, :], -1] = fluxunit
        
    # PSF
    if gdatmodi.thisindxprop == gdat.indxproppsfp:
        # temp
        if gdat.evalpsfnpnts:
            if gdat.modlvarioaxi:
                for k in gdat.indxoaxi:
                    gdatmodi.thispsfnintp[k] = interp1d(gdat.binsangl, gdatmodi.nextpsfn[:, :, :, k], axis=1)
            else:
                gdatmodi.thispsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
            
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
     
    # proposals that change the PS flux map
    if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfp:
        gdatmodi.thispntsflux[gdatmodi.indxcubemodi] = copy(gdatmodi.nextpntsflux[gdatmodi.indxcubemodi])

    # transdimensinal updates
    if gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]
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
        if gdat.numbener > 1:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompspep[gdatmodi.indxpoplmodi]]] = gdatmodi.modispep[0, gdat.indxspep[gdatmodi.indxpoplmodi]]
        
    ## death
    if gdatmodi.thisindxprop == gdat.indxpropdeth:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdatmodi.dethindxpnts)
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdatmodi.dethindxpnts)

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
        gdatmodi.thissampvarb[gdatmodi.indxsampfrst+2+gdat.numbener] = gdatmodi.modispep[1, gdat.indxspep[gdatmodi.indxpoplmodi]]
        #### second component
        gdatmodi.thissampvarb[gdatmodi.indxsampseco] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampseco+1] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampseco+2:gdatmodi.indxsampseco+2+gdat.numbener] = gdatmodi.modispec[:, 2]
        gdatmodi.thissampvarb[gdatmodi.indxsampseco+2+gdat.numbener] = gdatmodi.modispep[2, gdat.indxspep[gdatmodi.indxpoplmodi]]
        
    ## merge
    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdatmodi.mergindxseco)
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdatmodi.mergindxseco)

        ### update the component
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcomplgal]] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompbgal]] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompspec]] = gdatmodi.modispec[:, 2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompspep[gdatmodi.indxpoplmodi]]] = gdatmodi.modispep[2, gdat.indxspep[gdatmodi.indxpoplmodi]]
        
    ## PS parameter proposals

    if gdatmodi.thisindxprop >= gdat.indxproplgal:  
        if gdatmodi.thisindxprop == gdat.indxproplgal:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modilgal[1]
        elif gdatmodi.thisindxprop == gdat.indxpropbgal:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modibgal[1]
        else:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodispec] = gdatmodi.modispec[:, 1]
            if gdatmodi.thisindxprop == gdat.indxpropspep:
                if gdatmodi.indxspepmodi == gdat.indxspepsind:
                    gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modispep[1, 0]
                if gdatmodi.indxspepmodi == gdat.indxspepcurv or gdatmodi.indxspepmodi == gdat.indxspepexpo:
                    gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modispep[1, 1]


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


def retr_chandata(gdat):

    with open(gdat.pathdata + 'inpt/chancatl.txt', 'r') as thisfile:
        G_long = [] #deg
        G_lat = [] #deg
        id_number = []
        off_angle = [] # arcmin
        flux_cnts = [] # for xray band
        soft_cnts = [] # 0.5-2
        hard_cnts = [] # 2-8
        c_offset = [] #angular offset between xray and optical/NIR components in arcse
        C_mag = [] # full optical mag?
        W_mag = [] # 600-700 nm
        GD_mag = [] # 750-1000 nm from GOODS-S z-band
        G_mag = [] # 750-1000 nm from GEMS z-band
        M_mag = [] # 2-3 micron
        S_mag = [] # 3-4 micron
        flux_erg_full = [] # erg/(s*cm^2)
        flux_erg_soft = []
        flux_erg_hard = []
        Otype = [] # AGN/Galaxy/Star
        for line in thisfile:
            line = line.split()
            G_long.append(line[0])
            G_lat.append(line[1])
            id_number.append(line[2])
            off_angle.append(line[3])
            flux_cnts.append(line[4])
            soft_cnts.append(line[5])
            hard_cnts.append(line[6])
            c_offset.append(line[7])
            C_mag.append(line[8])
            W_mag.append(line[9])
            GD_mag.append(line[10])
            G_mag.append(line[11])
            M_mag.append(line[12])
            S_mag.append(line[13])
            flux_erg_full.append(line[14])
            flux_erg_soft.append(line[15])
            flux_erg_hard.append(line[16])
            Otype.append(line[17])
        lgalchan = (asarray(G_long)).astype(float)
        bgalchan = (asarray(G_lat)).astype(float)
        #oaxichan = (asarray(off_angle)).astype(float)
        cntschan = (asarray(flux_cnts)).astype(float)
        cntschansoft = (asarray(soft_cnts)).astype(float)
        cntschanhard = (asarray(hard_cnts)).astype(float)
        #offschan = (asarray(c_offset)).astype(float)
        #cmagchan = (asarray(C_mag)).astype(float)
        #wmagchan = (asarray(W_mag)).astype(float)
        #dmagchan = (asarray(GD_mag)).astype(float)
        #gmagchan = (asarray(G_mag)).astype(float)
        #mmagchan = (asarray(M_mag)).astype(float)
        #smagchan = (asarray(S_mag)).astype(float)
        #fluxchanfull = (asarray(flux_erg_full)).astype(float)
        fluxchansoft = (asarray(flux_erg_soft)).astype(float)
        fluxchanhard = (asarray(flux_erg_hard)).astype(float)
        #objttypechan = (asarray(Otype))

    #path = '/Users/tansu/Documents/work/data/tdgu/xray_back/data/CDFS-4Ms-0p5to2-bin1.emap'
    path = '/Users/tansu/Documents/work/data/tdgu/xray_back/data/CDFS-4Ms-0p5to2-asca-im-bin1.fits'
    listhdun = ap.io.fits.open(path)
    wcso = ap.wcs.WCS(listhdun[0].header)
    wx, wy = wcso.wcs_pix2world(1200, 1000, 1)
    a, b = wcso.wcs_world2pix(53.20, -27.9, 1)

    print 'w.wcs_pix2world(0, 0, 1)'
    print wcso.wcs_pix2world(0, 0, 1)
    print 'w.wcs_pix2world(0, 3279, 1)'
    print wcso.wcs_pix2world(0, 3279, 1)
    print 'w.wcs_pix2world(3259, 0, 1)'
    print wcso.wcs_pix2world(3259, 0, 1)
    print 'w.wcs_pix2world(0, 3280, 1)'
    print wcso.wcs_pix2world(0, 3280, 1)
    a, b = wcso.wcs_world2pix(lgalchan, bgalchan, 1)
    print a
    print b


    gdat.exprlgal = deg2rad(lgalchan)
    gdat.exprbgal = deg2rad(bgalchan)
    
    gdat.exprspec = zeros((3, gdat.numbener, gdat.exprlgal.size))
    gdat.exprcnts = zeros((gdat.numbener, gdat.exprlgal.size, gdat.numbevtt))
    gdat.exprspep = None

    gdat.exprspec[0, 0, :] = fluxchansoft * 0.624e12
    gdat.exprspec[0, 1, :] = fluxchanhard * 0.624e12 / 16.
    
    # temp
    gdat.exprspec[where(gdat.exprspec < 0.)] = 0.

    gdat.exprcnts[0, :, 0] = cntschansoft
    gdat.exprcnts[1, :, 0] = cntschanhard
    
    #gdat.exprstrg = lgalstrg
    #gdat.exprstrgclss = lgalchanclss
    #gdat.exprstrgassc = lgalchanassc

    indxsort = argsort(fluxchansoft)[::-1]
    
    gdat.exprlgal = gdat.exprlgal[indxsort][:150]
    gdat.exprbgal = gdat.exprbgal[indxsort][:150]
    gdat.exprspec = gdat.exprspec[:, :, indxsort][:150]
    gdat.exprcnts = gdat.exprcnts[:, indxsort][:150]


def retr_fermdata(gdat):
    
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'
    fgl3 = pf.getdata(path)
   
    gdat.exprlgal = deg2rad(fgl3['glon'])
    gdat.exprlgal = ((gdat.exprlgal - pi) % (2. * pi)) - pi
    gdat.exprbgal = deg2rad(fgl3['glat'])
    
    gdat.exprspec = empty((3, gdat.numbener, gdat.exprlgal.size))
    gdat.exprspec[0, :, :] = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], \
                                                                                            fgl3['Flux10000_100000']))[gdat.indxenerincl, :] / gdat.diffener[:, None]
    
    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.diffener[:, None, None] 
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :] - fgl3specstdvtemp[:, :, 0]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.exprspec[where(isfinite(gdat.exprspec) == False)] = 0.

    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3strg = fgl3['Source_Name']
    fgl3strgclss = fgl3['CLASS1']
    fgl3strgassc = fgl3['ASSOC1']
    
    gdat.exprspep = zeros((gdat.exprlgal.size, gdat.numbspeptotl))
    fgl3spectype = fgl3['SpectrumType']
    gdat.exprspep[:, 0] = fgl3['Spectral_Index']
    gdat.exprspep[:, 1] = fgl3['beta']
    gdat.exprspep[:, 2] = fgl3['Cutoff'] * 1e-3
    
    path = gdat.pathimag + '3fgl/'
    os.system('mkdir -p %s' % path)
    path += '3fglspecaxisstdv.pdf' 
    if not os.path.isfile(path):
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

        
def retr_angldistunit(gdat, lgal, bgal, indxpixltemp, retranglcosi=False):
   
    if gdat.pixltype == 'heal':
        xaxi, yaxi, zaxi = retr_unit(lgal, bgal)
        anglcosi = gdat.xaxigrid[indxpixltemp] * xaxi + gdat.yaxigrid[indxpixltemp] * yaxi + gdat.zaxigrid[indxpixltemp] * zaxi
        
        if retranglcosi:
            return anglcosi
        else:
            angldist = arccos(anglcosi)
            return angldist
    
    else:
        angldist = sqrt((lgal - gdat.lgalgrid[indxpixltemp])**2 + (bgal - gdat.bgalgrid[indxpixltemp])**2)
        
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


def show_samp(gdat, gdatmodi):
    print 'drmcsamp, thissampvarb'
    for k in range(gdatmodi.thissampvarb.size):
        if k == gdat.indxsampcomp[0]:
            print
        print '%14.4g %14.4g %14.4g' % (gdatmodi.drmcsamp[k, 0], gdatmodi.drmcsamp[k, 1], gdatmodi.thissampvarb[k])
    print


def retr_prop(gdat, gdatmodi):
 
    if gdat.verbtype > 1:
        print 'retr_prop(): '

        show_samp(gdat, gdatmodi)

        print 'listindxpntsfull: ', gdatmodi.listindxpntsfull
        print 'thisindxpntsfull: ', gdatmodi.thisindxpntsfull
        print 'thisindxpntsempt: ', gdatmodi.thisindxpntsempt  
        print 'thisindxsamplgal: ', gdatmodi.thisindxsamplgal
        print 'thisindxsampbgal: ', gdatmodi.thisindxsampbgal
        print 'thisindxsampspec: '
        print gdatmodi.thisindxsampspec
        print 'thisindxsampspep: ', gdatmodi.thisindxsampspep
        print 'thisindxsampcompcolr: ', gdatmodi.thisindxsampcompcolr
        print
        
    # fixed dimensional parameter change 
    if gdatmodi.propfixp:
    
        # take the step
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvproppara[gdatmodi.indxsampmodi])

        ## hyperparameter changes
        if gdatmodi.thisindxprop == gdat.indxpropfluxdist or gdatmodi.thisindxprop == gdat.indxpropmeanpnts:
            gdatmodi.indxpoplmodi = (gdatmodi.indxsampmodi - gdat.indxfixpfluxdist[0]) // gdat.numbfluxdistpara
            
        # PSF parameter changes 
        if gdatmodi.thisindxprop == gdat.indxproppsfp:
            gdatmodi.nextsampvarb[gdat.indxfixppsfp] = copy(gdatmodi.thissampvarb[gdat.indxfixppsfp])
                
            ## determine the background index to be modified
            gdatmodi.indxenermodi = array([((gdatmodi.indxsampmodi - gdat.indxfixppsfp[0]) % gdat.numbpsfptotlener) // gdat.numbpsfptotl])
            
            if gdat.verbtype > 1:
                print 'thispsfp'
                print gdatmodi.thissampvarb[gdat.indxfixppsfp]
                print 'nextpsfp'
                print gdatmodi.nextsampvarb[gdat.indxfixppsfp]
        
        # background parameter changes
        if gdatmodi.thisindxprop == gdat.indxpropbacp:

            ## determine the background index to be modified
            gdatmodi.indxenermodi = (gdatmodi.indxsampmodi - gdat.indxfixpbacp[0]) // gdat.numbback
            gdatmodi.indxbackmodi = (gdatmodi.indxsampmodi - gdat.indxfixpbacp[0]) % gdat.numbback
            
            ## save the current background parameters
            gdatmodi.nextsampvarb[gdat.indxfixpbacp] = copy(gdatmodi.thissampvarb[gdat.indxfixpbacp])
            
            if gdat.verbtype > 1:
                print 'indxbackmodi'
                print gdatmodi.indxbackmodi
                print 'indxenermodi'
                print gdatmodi.indxenermodi

        # lens parameter changes
        if gdatmodi.thisindxprop == gdat.indxproplenp:
            # temp
            gdatmodi.indxenermodi = 0

        # inverse CDF transform the proposed variable 
        gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_fixp(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdatmodi.indxsampmodi)
    
    # birth
    if gdatmodi.thisindxprop == gdat.indxpropbrth:

        # temp -- modi
        gdatmodi.numbpntsmodi = 1
        #thisnumbpntsmodi = gdat.maxmnumbpnts[gdatmodi.indxpoplmodi] - int(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]])
        #gdatmodi.numbpntsmodi = choice(gdat.listnumbpntsmodi[thisnumbpntsmodi], p=gdat.probnumbpntsmodi[thisnumbpntsmodi])

        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + gdatmodi.numbpntsmodi
    
        # initial sample index to add the new PS
        # temp -- modi
        indxsampbrth = gdat.indxsampcomp[0] + gdat.maxmnumbcompcuml[gdatmodi.indxpoplmodi] + \
                                                       array(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][:gdatmodi.numbpntsmodi]) * gdat.numbcomp[gdatmodi.indxpoplmodi]
        
        # sample auxiliary variables
        # number of sample vector elements to be modified
        numbcompmodi = gdatmodi.numbpntsmodi * gdat.numbcomp[gdatmodi.indxpoplmodi]
        # number of unit sample vector elements to be modified
        numbcompcolrmodi = gdatmodi.numbpntsmodi * gdat.numbcompcolr[gdatmodi.indxpoplmodi]
        # auxiliary vector
        gdatmodi.auxipara = rand(numbcompcolrmodi)
        # index of samples to be modified
        gdatmodi.indxsampmodi = empty(numbcompmodi, dtype=int)
        for k in range(gdatmodi.numbpntsmodi):
            gdatmodi.drmcsamp[indxsampbrth[k]+gdat.indxcompcolr[gdatmodi.indxpoplmodi], -1] = gdatmodi.auxipara[k*gdatmodi.numbpntsmodi+gdat.indxauxipara[gdatmodi.indxpoplmodi]]
            # sample indices to be modified
            gdatmodi.indxsampmodi[k*gdat.numbcomp[gdatmodi.indxpoplmodi]:(k+1)*gdat.numbcomp[gdatmodi.indxpoplmodi]] = indxsampbrth[k] + gdat.indxcomp[gdatmodi.indxpoplmodi]

        # modification catalog
        gdatmodi.modilgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcomplgal, -1], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
        gdatmodi.modibgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompbgal, -1], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
        fluxunit = gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompflux, -1]
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
            fluxdistslop = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]
            gdatmodi.modispec[gdat.indxenerfluxdist, :gdatmodi.numbpntsmodi] = icdf_flux_powr(fluxunit, gdat.minmflux, gdat.maxmflux, fluxdistslop)
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
            fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]]
            fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]
            fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]]
            gdatmodi.modispec[gdat.indxenerfluxdist, :gdatmodi.numbpntsmodi] = icdf_flux_brok(array([fluxunit]), gdat.minmflux, gdat.maxmflux, \
                                                                                                        fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        if gdat.numbener > 1:
            gdatmodi.modispep[0, gdat.indxspep[gdatmodi.indxpoplmodi]] = icdf_gaus(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompsind, -1], \
                                                                                gdat.sinddistmean[gdatmodi.indxpoplmodi], gdat.sinddiststdv[gdatmodi.indxpoplmodi])
            gdatmodi.modispec[:, 0] = retr_spec(gdat, gdatmodi.modispec[gdat.indxenerfluxdist, 0], \
                                                                                                    spep=gdatmodi.modispep[0, gdat.indxspep[gdatmodi.indxpoplmodi]], \
                                                                                                    spectype=gdat.spectype[gdatmodi.indxpoplmodi]).flatten()
    
        if gdat.verbtype > 1:
            print 'numbpntsmodi'
            print gdatmodi.numbpntsmodi
            print 'auxipara'
            print gdatmodi.auxipara
            print 'numbcompcolrmodi'
            print numbcompcolrmodi
            print 'numbcompmodi'
            print numbcompmodi
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
    # death
    if gdatmodi.thisindxprop == gdat.indxpropdeth:
        
        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # occupied PS index to be killed
        dethindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        gdatmodi.dethindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][dethindxindxpnts]
        
        # sample indices to be modified 
        gdatmodi.indxsampmodi = array([])
            
        # modification catalog
        gdatmodi.numbpntsmodi = 1
        gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][dethindxindxpnts]]
        gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][dethindxindxpnts]]
        gdatmodi.modispec[:, 0] = -gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, dethindxindxpnts]]

        if gdat.verbtype > 1:
            print 'dethindxpnts: ', gdatmodi.dethindxpnts
            print 'dethindxindxpnts: ', dethindxindxpnts
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
  
    # split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:
        
        gdatmodi.numbpntsmodi = 3
        
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        gdatmodi.spltindxindxpnts = choice(thisindxindxpnts)
    
        # update the sample vector
        gdatmodi.indxsampfrst = gdat.indxsampcomp[0] + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + \
                                                int(gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]) * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp[gdatmodi.indxpoplmodi]

        print 'gdatmodi.thisindxpntsempt'
        print gdatmodi.thisindxpntsempt
        print 'gdatmodi.indxpoplmodi'
        print gdatmodi.indxpoplmodi
        print 'gdat.indxsampcomp[0]'
        print gdat.indxsampcomp[0]
        print 
        gdatmodi.indxsampseco = gdat.indxsampcomp[0] + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + \
                                                int(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0]) * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlseco = gdatmodi.indxsampseco + gdat.numbcomp[gdatmodi.indxpoplmodi]
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, indxfinlfrst, dtype=int), arange(gdatmodi.indxsampseco, indxfinlseco, dtype=int)))
        
        thislgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisbgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, gdatmodi.spltindxindxpnts]]
        gdatmodi.fluxpare = thisspec[gdat.indxenerfluxdist[0]]
        thisspep = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts, :]]
        
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
            print 'thisspep: ', thisspep
            
        # determine the new components
        # temp -- only valid for power-law energy spectrum
        gdatmodi.auxipara = empty(gdat.numbcompcolr[gdatmodi.indxpoplmodi])
        gdatmodi.auxipara[0] = rand()
        gdatmodi.auxipara[1] = rand() * gdat.radispmr
        gdatmodi.auxipara[2] = rand() * 2. * pi
        # temp
        if gdat.numbener > 1:
            gdatmodi.auxipara[3] = icdf_gaus(rand(), gdat.sinddistmean[gdatmodi.indxpoplmodi], gdat.sinddiststdv[gdatmodi.indxpoplmodi])
            

        if gdat.verbtype > 1:
            print 'auxipara[0]: ', gdatmodi.auxipara[0]
            print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
            print 'auxipara[2]: ', gdatmodi.auxipara[2]
            if gdat.numbener > 1:
                print 'auxipara[3]: ', gdatmodi.auxipara[3]
            print
            
        gdatmodi.fluxfrst = gdatmodi.auxipara[0] * gdatmodi.fluxpare
        gdatmodi.spltlgalfrst = thislgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1] * cos(gdatmodi.auxipara[2])
        gdatmodi.spltbgalfrst = thisbgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1] * sin(gdatmodi.auxipara[2])
        gdatmodi.spltsindfrst = thisspep
        
        gdatmodi.fluxseco = (1. - gdatmodi.auxipara[0]) * gdatmodi.fluxpare
        gdatmodi.spltlgalseco = thislgal - gdatmodi.auxipara[0] * gdatmodi.auxipara[1] * cos(gdatmodi.auxipara[2])
        gdatmodi.spltbgalseco = thisbgal - gdatmodi.auxipara[0] * gdatmodi.auxipara[1] * sin(gdatmodi.auxipara[2])
        if gdat.numbener > 1:
            gdatmodi.spltsindseco = gdatmodi.auxipara[3]
        
        if gdat.verbtype > 1:
            print 'spltlgalfrst: ', gdat.anglfact * gdatmodi.spltlgalfrst
            print 'spltlgalseco: ', gdat.anglfact * gdatmodi.spltlgalseco
            print 'spltbgalfrst: ', gdat.anglfact * gdatmodi.spltbgalfrst
            print 'spltbgalseco: ', gdat.anglfact * gdatmodi.spltbgalseco
            print 'spltfluxfrst: ', gdatmodi.fluxfrst
            print 'spltfluxseco: ', gdatmodi.fluxseco
            if gdat.numbener > 1:
                print 'spltsindfrst: ', gdatmodi.spltsindfrst
                print 'spltsindseco: ', gdatmodi.spltsindseco

        if fabs(gdatmodi.spltlgalfrst) > gdat.maxmgangmodl or fabs(gdatmodi.spltlgalseco) > gdat.maxmgangmodl or \
                                            fabs(gdatmodi.spltbgalfrst) > gdat.maxmgangmodl or fabs(gdatmodi.spltbgalseco) > gdat.maxmgangmodl or \
                                            gdatmodi.fluxfrst < gdat.minmflux or gdatmodi.fluxseco < gdat.minmflux:
            gdatmodi.boolreje = True

        if gdat.verbtype > 1:
            print 'boolreje'
            print gdatmodi.boolreje

        # calculate the list of pairs
        ## current
        gdatmodi.thislistpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], \
                                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]])
        gdatmodi.thisnumbpair = len(gdatmodi.thislistpair)
        
        if not gdatmodi.boolreje:

            # calculate the list of pairs
            ## proposed
            lgal = concatenate((array([gdatmodi.spltlgalfrst, gdatmodi.spltlgalseco]), \
                                                                    setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([gdatmodi.spltbgalfrst, gdatmodi.spltbgalseco]), \
                                                                    setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]], thisbgal)))
            gdatmodi.nextlistpair = retr_listpair(gdat, lgal, bgal)
            gdatmodi.nextnumbpair = len(gdatmodi.nextlistpair)

            if gdatmodi.nextnumbpair == 0:
                print 'Number of pairs should not be zero in the reverse proposal of a split'
                raise

            ## first new component
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcomplgal, -1] = cdfn_self(gdatmodi.spltlgalfrst, -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcompbgal, -1] = cdfn_self(gdatmodi.spltbgalfrst, -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcompflux, -1] = cdfn_flux_powr(gdatmodi.fluxfrst, gdat.minmflux, gdat.maxmflux, \
                                                                                    gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+gdat.indxcompsind, -1] = cdfn_gaus(gdatmodi.spltsindfrst, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
    
            # make retr_spec be called only for gdat.numbener > 1
            nextspecfrst = retr_spec(gdat, gdatmodi.fluxfrst, spep=gdatmodi.spltsindfrst, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

            ## second new component
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcomplgal, -1] = cdfn_self(gdatmodi.spltlgalseco, -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcompbgal, -1] = cdfn_self(gdatmodi.spltbgalseco, -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcompflux, -1] = cdfn_flux_powr(gdatmodi.fluxseco, gdat.minmflux, gdat.maxmflux, \
                                                                                    gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdatmodi.indxsampseco+gdat.indxcompsind, -1] = cdfn_gaus(gdatmodi.spltsindseco, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                                                                                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
            nextspecseco = retr_spec(gdat, gdatmodi.fluxseco, spep=gdatmodi.spltsindseco, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

            ## component to be removed
            gdatmodi.modilgal[0] = thislgal
            gdatmodi.modibgal[0] = thisbgal
            gdatmodi.modispec[:, 0] = -thisspec.flatten()
            gdatmodi.modispep[0, gdat.indxspep[gdatmodi.indxpoplmodi]] = thisspep
            
            ## first component to be added
            gdatmodi.modilgal[1] = gdatmodi.spltlgalfrst
            gdatmodi.modibgal[1] = gdatmodi.spltbgalfrst
            gdatmodi.modispec[:, 1] = nextspecfrst.flatten()
            gdatmodi.modispep[1, gdat.indxspep[gdatmodi.indxpoplmodi]] = gdatmodi.spltsindfrst

            # second component to be added
            gdatmodi.modilgal[2] = gdatmodi.spltlgalseco
            gdatmodi.modibgal[2] = gdatmodi.spltbgalseco
            gdatmodi.modispec[:, 2] = nextspecseco.flatten()
            gdatmodi.modispep[2, gdat.indxspep[gdatmodi.indxpoplmodi]] = gdatmodi.spltsindseco

    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        # number of point sources to be modified
        gdatmodi.numbpntsmodi = 3
        
        # proposed number of point sources
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # calculate the current list of pairs
        gdatmodi.thislistpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]])
        gdatmodi.thisnumbpair = len(gdatmodi.thislistpair)
        if gdat.verbtype > 1:
            print 'thislistpair'
            print gdatmodi.thislistpair
           
        # check if merge will be proposed
        if gdatmodi.thisnumbpair == 0:
            gdatmodi.boolreje = True
        else:

            # sample a pair
            indxpairtemp = choice(arange(gdatmodi.thisnumbpair))

            # determine PS indices to be merged
            mergindxindxpntsfrst = gdatmodi.thislistpair[indxpairtemp][0]
            mergindxindxpntsseco = gdatmodi.thislistpair[indxpairtemp][1]
  
            ## first PS index to be merged
            gdatmodi.mergindxfrst = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]

            ## second PS index to be merged
            gdatmodi.mergindxseco = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpntsseco]

            # determine indices of the modified elements in the sample vector
            ## first PS
            # temp -- this would not work for multiple populations !
            gdatmodi.indxsampfrst = gdat.indxsampcomp[0] + gdat.numbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxfrst
            indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp[gdatmodi.indxpoplmodi]
            
            ## second PS
            gdatmodi.indxsampseco = gdat.indxsampcomp[0] + gdat.numbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxseco
            indxfinlseco = gdatmodi.indxsampseco + gdat.numbcomp[gdatmodi.indxpoplmodi]

            # indices of the sample vector elements to be modified
            gdatmodi.indxsampmodi = arange(gdatmodi.indxsampfrst, indxfinlfrst)

            # indices of the PS to be merged
            mergindxpnts = sort(array([gdatmodi.mergindxfrst, gdatmodi.mergindxseco], dtype=int))

            # PS parameters to be merged
            ## first PS
            gdatmodi.lgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.bgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.specfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpntsfrst]]
            gdatmodi.spepfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][mergindxindxpntsfrst, :]]
            gdatmodi.fluxfrst = gdatmodi.specfrst[gdat.indxenerfluxdist[0]]

            ## second PS
            gdatmodi.lgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.bgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.specseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpntsseco]]
            gdatmodi.spepseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][mergindxindxpntsseco, :]]
            gdatmodi.fluxseco = gdatmodi.specseco[gdat.indxenerfluxdist[0]]

            # auxiliary parameters
            auxifrac = gdatmodi.fluxfrst / (gdatmodi.fluxfrst + gdatmodi.fluxseco) 
            auxiradi = sqrt((gdatmodi.lgalseco - gdatmodi.lgalfrst)**2 + (gdatmodi.bgalseco - gdatmodi.bgalfrst)**2)
            auxiangl = pi + arctan2(gdatmodi.bgalseco - gdatmodi.bgalfrst, gdatmodi.lgalseco - gdatmodi.lgalfrst)
            auxispep = gdatmodi.spepseco

            # temp
            gdatmodi.auxipara = zeros(gdat.numbcompcolr[gdatmodi.indxpoplmodi])
            gdatmodi.auxipara[0] = auxifrac
            gdatmodi.auxipara[1] = auxiradi
            gdatmodi.auxipara[2] = auxiangl
            gdatmodi.auxipara[3:] = gdatmodi.spepseco
            
            # merged PS
            gdatmodi.fluxpare = gdatmodi.fluxfrst + gdatmodi.fluxseco
            if gdatmodi.fluxpare > gdat.maxmflux:
                gdatmodi.boolreje = True
            gdatmodi.lgalpare = gdatmodi.lgalfrst + (1. - auxifrac) * (gdatmodi.lgalseco - gdatmodi.lgalfrst)
            gdatmodi.bgalpare = gdatmodi.bgalfrst + (1. - auxifrac) * (gdatmodi.bgalseco - gdatmodi.bgalfrst)
            gdatmodi.speppare = gdatmodi.spepfrst
            gdatmodi.specpare = retr_spec(gdat, gdatmodi.fluxpare, spep=gdatmodi.speppare, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

            # determine the unit variables for the merged PS
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst, -1] = cdfn_self(gdatmodi.lgalpare, -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+1, -1] = cdfn_self(gdatmodi.bgalpare, -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+2, -1] = cdfn_flux_powr(gdatmodi.fluxpare, gdat.minmflux, gdat.maxmflux, \
                                                                                            gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdatmodi.indxsampfrst+3, -1] = gdatmodi.drmcsamp[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][mergindxindxpntsfrst, :], -2]

            # PSs to be added to the PS flux map
            ## first PS
            gdatmodi.modilgal[0] = gdatmodi.lgalfrst
            gdatmodi.modibgal[0] = gdatmodi.bgalfrst
            gdatmodi.modispec[:, 0] = -gdatmodi.specfrst.flatten()
            gdatmodi.modispep[0, gdat.indxspep[gdatmodi.indxpoplmodi]] = gdatmodi.spepfrst

            ## second PS
            gdatmodi.modilgal[1] = gdatmodi.lgalseco
            gdatmodi.modibgal[1] = gdatmodi.bgalseco
            gdatmodi.modispec[:, 1] = -gdatmodi.specseco.flatten()
            gdatmodi.modispep[1, gdat.indxspep[gdatmodi.indxpoplmodi]] = gdatmodi.spepseco

            ## parent PS
            gdatmodi.modilgal[2] = gdatmodi.lgalpare
            gdatmodi.modibgal[2] = gdatmodi.bgalpare
            gdatmodi.modispec[:, 2] = gdatmodi.specpare.flatten()
            gdatmodi.modispep[2, gdat.indxspep[gdatmodi.indxpoplmodi]] = gdatmodi.speppare

            # calculate the proposed list of pairs
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
                print 'mergspepfrst: ', gdatmodi.spepfrst
                print 'merglgalseco: ', gdat.anglfact * gdatmodi.lgalseco
                print 'mergbgalseco: ', gdat.anglfact * gdatmodi.bgalseco
                print 'mergfluxseco: ', gdatmodi.fluxseco
                print 'mergspepseco: ', gdatmodi.spepseco
                print 'merglgalpare: ', gdat.anglfact * gdatmodi.lgalpare
                print 'mergbgalpare: ', gdat.anglfact * gdatmodi.bgalpare
                print 'mergspecpare: ', gdatmodi.specpare
                print 'mergfluxpare: ', gdatmodi.fluxpare
                print 'mergspeppare: ', gdatmodi.speppare
                print 'auxipara[0]: ', gdatmodi.auxipara[0]
                print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
                print 'auxipara[2]: ', gdatmodi.auxipara[2]
                if gdat.numbener > 1:
                    print 'auxipara[3]: ', gdatmodi.auxipara[3]
            
            lgal = concatenate((array([gdatmodi.lgalpare]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], \
                                                                                                        array([gdatmodi.lgalfrst, gdatmodi.lgalseco]))))
            bgal = concatenate((array([gdatmodi.bgalpare]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]], \
                                                                                                        array([gdatmodi.bgalfrst, gdatmodi.bgalseco]))))
            gdatmodi.nextlistpair = retr_listpair(gdat, lgal, bgal)
            gdatmodi.nextnumbpair = len(gdatmodi.nextlistpair)
        
            if gdat.verbtype > 1:
                # temp
                if False:
                    print 'nextlistpair'
                    print gdatmodi.nextlistpair
                print

            if auxiradi > gdat.radispmr:
                print 'Auxiliary radius during a merge cannot be larger than the linking length of %.3g %s.' % (gdat.anglfact * gdat.radispmr, gdat.strganglunit)
                raise

    # PS parameter change
    if gdatmodi.thisindxprop >= gdat.indxproplgal:     
        
        gdatmodi.indxenermodi = gdat.indxener
            
        # initial sample index of the PS to be modified
        gdatmodi.indxsampmodiinit = gdat.indxsampcomp[0] + gdat.maxmnumbcompcuml[gdatmodi.indxpoplmodi] + gdatmodi.indxpntsmodi * gdat.numbcomp[gdatmodi.indxpoplmodi]
        
        # sample index to be modified
        gdatmodi.indxsampmodispec = gdatmodi.indxsampmodiinit + 2 + gdat.indxener
        
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, gdatmodi.indxpntsfullmodi]]
        if gdat.numbener > 1:
            thisspep = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, :]]
            thisspec = retr_spec(gdat, thisspec[gdat.indxenerfluxdist], spep=thisspep, spectype=gdat.spectype[gdatmodi.indxpoplmodi])
            
        # propose
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdat.stdvlbhlvari:
                retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvlbhl * gdat.minmflux / thisspec[gdat.indxenerfluxdist])
            else:
                retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvlbhl) 
        if gdatmodi.thisindxprop == gdat.indxpropflux:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvflux)
        if gdatmodi.thisindxprop == gdat.indxpropspep:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvspep)

        gdatmodi.numbpntsmodi = 2
        gdatmodi.modispec[:, 0] = -thisspec.flatten()
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdatmodi.indxcompmodi == 0:
                gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi]]
                gdatmodi.modilgal[1] = icdf_self(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
                gdatmodi.modibgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi]]
            else:
                gdatmodi.modilgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi]]
                gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi]]
                gdatmodi.modibgal[1] = icdf_self(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            gdatmodi.modispec[:, 1] = thisspec.flatten()
        else:
            gdatmodi.modilgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi]]
            gdatmodi.modibgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi]]
            if gdatmodi.thisindxprop == gdat.indxpropflux:
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]
                    gdatmodi.modispec[gdat.indxenerfluxdist, :gdatmodi.numbpntsmodi] = icdf_flux_powr(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
                                                                                                                        gdat.minmflux, gdat.maxmflux, fluxdistslop)
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                    gdatmodi.modispec[gdat.indxenerfluxdist, :gdatmodi.numbpntsmodi] = icdf_flux_brok(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.minmflux, \
                                                                                                        gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
                
                if gdat.numbener > 1:
                    gdatmodi.modispep[1, 0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, 0]]
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'curv':
                        gdatmodi.modispep[1, 1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, 1]]
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'expo':
                        gdatmodi.modispep[1, 1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, 1]]
            else:
                gdatmodi.modispec[gdat.indxenerfluxdist, :gdatmodi.numbpntsmodi] = \
                                            gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, gdatmodi.indxpntsfullmodi]]
                if gdatmodi.indxspepmodi == gdat.indxspepsind:
                    gdatmodi.modispep[1, 0] = icdf_gaus(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                                                                                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'curv':
                        gdatmodi.modicurv[1, 1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, 1]]
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'expo':
                        gdatmodi.modispep[1, 1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, 1]]
           
                else:
                    gdatmodi.modispep[1, 0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.indxpntsfullmodi, 0]]
                    if gdatmodi.indxspepmodi == gdat.indxspepcurv:
                        gdatmodi.modispep[1, 1] = icdf_gaus(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                                                                                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
                    if gdatmodi.indxspepmodi == gdat.indxspepexpo:
                        gdatmodi.modispep[1, 1] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.minmener, gdat.factener)
           
            if gdat.numbener > 1:
                gdatmodi.modispec[:, 1] = retr_spec(gdat, gdatmodi.modispec[gdat.indxenerfluxdist[0], 1], \
                                                                                                        spep=gdatmodi.modispep[1, gdat.indxspep[gdatmodi.indxpoplmodi]], \
                                                                                                        spectype=gdat.spectype[gdatmodi.indxpoplmodi]).flatten()

        if gdat.verbtype > 1:
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print 'indxcompmodi: ', gdatmodi.indxcompmodi
            print 'indxpntsfullmodi: ', gdatmodi.indxpntsfullmodi
            print 'indxpntsmodi'
            print gdatmodi.indxpntsmodi

    # energy bin in which to evaluate the log-likelihood
    if gdat.indxpropbrth <= gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdatmodi.indxenermodi = gdat.indxener

    if gdat.verbtype > 1:
        if gdat.numbener > 1 and gdatmodi.thisindxprop in gdat.indxpropllik:
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
        thisnumbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]
        if gdatmodi.thisindxprop == gdat.indxpropsplt:
            gdatmodi.combfact = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2 / gdatmodi.nextnumbpair
        else:
            gdatmodi.combfact = gdatmodi.thisnumbpair / gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2
        
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
   
        
def retr_factoaxi(gdat, bins, norm, indx):

    factoaxi = 1. + norm[:, None] * (bins[None, :] / gdat.oaxipivt)**indx[:, None]
     
    return factoaxi


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, binsoaxi=None, varioaxi=None):

    numbpsfpform, numbpsfpoaxi = retr_numbpsfp(psfntype, varioaxi)
    print 'numbpsfpform'
    print numbpsfpform
    print 'numbpsfpoaxi'
    print numbpsfpoaxi
    numbpsfptotl = numbpsfpform + numbpsfpoaxi
    indxpsfptemp = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    print 'indxpsfptemp'
    print indxpsfptemp
    if varioaxi:
        indxpsfpoaxinorm = numbpsfptotl * gdat.indxener[indxenertemp] + gdat.numbpsfpform
        indxpsfpoaxiindx = numbpsfptotl * gdat.indxener[indxenertemp] + gdat.numbpsfpform + 1
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * arcsin(sqrt(2. - 2. * cos(gdat.binsangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        if varioaxi:
            scalangl = thisangl[None, :, None, None]
        else:
            scalangl = thisangl[None, :, None]
    
    if varioaxi:
        factoaxi = retr_factoaxi(gdat, binsoaxi, psfp[indxpsfpoaxinorm], psfp[indxpsfpoaxiindx])
   
        print 'factoaxi'
        print factoaxi
        print 

    if psfntype == 'singgaus':
        sigc = psfp[indxpsfptemp]
        if varioaxi:
            sigc = sigc[:, None, :, None] * factoaxi[:, None, None, :]
        else:
            sigc = sigc[:, None, :]

        print 'sigc'
        print sigc[0, :, 0, 0]
        print sigc[1, :, 0, 0]
        print 'scalangl'
        print scalangl[0, :, 0, 0]
        print 

        psfn = retr_singgaus(scalangl, sigc)
        
    elif psfntype == 'singking':
        sigc = psfp[indxpsfptemp]
        gamc = psfp[indxpsfptemp+1]
        psfn = retr_singking(scalangl, sigc, gamc)
        if varioaxi:
            sigc = sigc[:, None, :, None] * factoaxi[None, None, None, :]
            gamc = gamc[:, None, :, None]
        else:
            sigc = sigc[:, None, :]
            gamc = gamc[:, None, :]
        
    elif psfntype == 'doubgaus':
        frac = psfp[indxpsfptemp]
        sigc = psfp[indxpsfptemp+1]
        sigt = psfp[indxpsfptemp+2]
        if varioaxi:
            frac = frac[:, None, :, None]
            sigc = sigc[:, None, :, None] * factoaxi[None, None, None, :]
            sigt = sigt[:, None, :, None] * factoaxi[None, None, None, :]
        else:
            frac = frac[:, None, :]
            sigc = sigc[:, None, :]
            sigt = sigt[:, None, :]
        psfn = retr_doubgaus(scalangl, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfp[indxpsfptemp]
        sigc = psfp[indxpsfptemp+1]
        sigt = psfp[indxpsfptemp+2]
        gamt = psfp[indxpsfptemp+3]
        if varioaxi:
            frac = frac[:, None, :, None]
            sigc = sigc[:, None, :, None] * factoaxi[None, None, None, :]
            sigt = sigt[:, None, :, None] * factoaxi[None, None, None, :]
            gamt = gamt[:, None, :, None]
        else:
            frac = frac[:, None, :]
            sigc = sigc[:, None, :]
            sigt = sigt[:, None, :]
            gamt = gamt[:, None, :]
        psfn = retr_gausking(scalangl, frac, sigc, sigt, gamt)
        
    elif psfntype == 'doubking':
        frac = psfp[indxpsfptemp]
        sigc = psfp[indxpsfptemp+1]
        gamc = psfp[indxpsfptemp+2]
        sigt = psfp[indxpsfptemp+3]
        gamt = psfp[indxpsfptemp+4]
        if varioaxi:
            frac = frac[:, None, :, None]
            sigc = sigc[:, None, :, None] * factoaxi[None, None, None, :]
            gamc = gamc[:, None, :, None]
            sigt = sigt[:, None, :, None] * factoaxi[None, None, None, :]
            gamt = gamt[:, None, :, None]
        else:
            frac = frac[:, None, :]
            sigc = sigc[:, None, :]
            gamc = gamc[:, None, :]
            sigt = sigt[:, None, :]
            gamt = gamt[:, None, :]

        psfn = retr_doubking(scalangl, frac, sigc, gamc, sigt, gamt)
        if gdat.exprtype == 'ferm':
            psfnnorm = retr_doubking(scalanglnorm, frac, sigc, gamc, sigt, gamt)
    
    # normalize the PSF
    if gdat.exprtype == 'ferm':
        fact = 2. * pi * trapz(psfnnorm * sin(gdat.binsangl[None, :, None]), gdat.binsangl, axis=1)[:, None, :]
        psfn /= fact

    # temp
    if True and (gdat.strgcnfg == 'pcat_ferm_expr_ngal' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp1' or \
                                                            gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp2' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp3'):
        print 'CORRECTING THE PSF.'
        tempcorr = array([1., 0.8, 0.8])
        psfn *= tempcorr[:, None, None]

    return psfn


def retr_unit(lgal, bgal):

    xaxi = cos(bgal) * cos(lgal)
    yaxi = -cos(bgal) * sin(lgal)
    zaxi = sin(bgal)

    return xaxi, yaxi, zaxi


def retr_randunitpsfp(gdat):

    while True:
        randunitpsfp = rand(gdat.numbpsfp)
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
                    indx = m * gdat.numbpsfpevtt + i * gdat.numbpsfpform
                    thisbool = thisbool and randunitpsfp[indx+indxpar1] > randunitpsfp[indx+indxpar0]
            if thisbool:
                break

    return randunitpsfp


def retr_eerrnorm(minmvarb, maxmvarb, meanvarb, stdvvarb):
   
    cdfnminm = 0.5 * (sp.special.erf((minmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfnmaxm = 0.5 * (sp.special.erf((maxmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfndiff = cdfnmaxm - cdfnminm
    
    return cdfnminm, cdfndiff


def retr_numbpsfp(psfntype, varioaxi):

    if psfntype == 'singgaus':
        numbpsfpform = 1
    elif psfntype == 'singking':
        numbpsfpform = 2
    elif psfntype == 'doubgaus':
        numbpsfpform = 3
    elif psfntype == 'gausking':
        numbpsfpform = 4
    elif psfntype == 'doubking':
        numbpsfpform = 5
    
    if varioaxi:
        numbpsfpoaxi = 2
    else:
        numbpsfpoaxi = 0

    return numbpsfpform, numbpsfpoaxi


def retr_numbspep(spectype):
    
    numbpopl = len(spectype)
    numbspep = empty(numbpopl, dtype=int)
    for l in range(numbpopl):
        if spectype[l] == 'powr':
            numbspep[l] = 1
        if spectype[l] == 'expo':
            numbspep[l] = 2
        if spectype[l] == 'curv':
            numbspep[l] = 2
     
    return numbspep
    

def setpinit(gdat, boolinitsetp=False):

    # paths
    gdat.pathdata = gdat.pathbase + 'data/'
    gdat.pathimag = gdat.pathbase + 'imag/'
    
    # process index
    gdat.indxproc = arange(gdat.numbproc)

    # population index vector
    if gdat.datatype == 'mock':
        gdat.mockindxpopl = arange(gdat.mocknumbpopl, dtype=int)

    # flag to indicate whether information from a deterministic catalog will be used or not
    # temp -- if datatype == 'inpt' trueinfo should depend on whether truexxxx are provided
    gdat.trueinfo = gdat.datatype == 'mock' or gdat.exprinfo
    
    # half size of the image where the sample catalog is compared against the reference
    gdat.maxmgangcomp = gdat.maxmgang * gdat.margfactcomp
    # half size of the spatial prior
    gdat.maxmgangmodl = gdat.maxmgang * gdat.margfactmodl

    # axes
    gdat.numbfluxplot = 20
    
    gdat.indxspepsind = 0
    gdat.indxspepcurv = 1
    gdat.indxspepexpo = 2

    gdat.numbfluxdistpara = 4
    
    if gdat.strgfluxunit == None:
        gdat.strgfluxunitextn = ''
    else:
        gdat.strgfluxunitextn = ' [%s]' % gdat.strgfluxunit

    ## energy
    if gdat.binsenerfull != None:
        gdat.enerbins = True
    else:
        gdat.enerbins = False
    
    if gdat.indxevttincl != None:
        gdat.evttbins = True
    else:
        gdat.evttbins = False

    if gdat.enerbins:
        gdat.numbener = gdat.indxenerincl.size
        gdat.numbenerfull = gdat.binsenerfull.size - 1
        gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
        gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
        gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
        gdat.binsener = gdat.binsenerfull[gdat.indxenerinclbins]
        gdat.diffener = (roll(gdat.binsener, -1) - gdat.binsener)[0:-1]
        gdat.meanener = sqrt(roll(gdat.binsener, -1) * gdat.binsener)[0:-1]
        gdat.minmener = gdat.binsener[0]
        gdat.maxmener = gdat.binsener[-1]
    else:
        gdat.numbener = 1
        gdat.numbenerfull = 1
        gdat.indxenerincl = array([0])
        gdat.indxenerfluxdist = array([0])
        gdat.factspecener = array([1.])
    gdat.indxener = arange(gdat.numbener, dtype=int)
    gdat.indxenerfluxdist = ceil(array([gdat.numbener]) / 2.).astype(int) - 1
       
    if gdat.enerbins:
        gdat.enerfluxdist = gdat.meanener[gdat.indxenerfluxdist]
        if gdat.enerfluxdist == 0.:
            raise Exception('Pivot energy cannot be zero.')
        gdat.enernorm = gdat.meanener / gdat.enerfluxdist
        gdat.factlogtenerpivt = log(gdat.enernorm)
        gdat.factspecener = gdat.enernorm**(-gdat.sinddistmean)
        gdat.enerexpofact = gdat.enerfluxdist - gdat.meanener

    ## PSF class
    if gdat.evttbins:
        gdat.numbevtt = gdat.indxevttincl.size
        gdat.numbevttfull = gdat.indxevttfull.size
    else:
        gdat.numbevtt = 1
        gdat.numbevttfull = 1
        gdat.indxevttincl = array([0])
    gdat.indxevtt = arange(gdat.numbevtt)

    # angular deviation
    # temp -- check that gdat.numbangl does not degrade the performance
    if gdat.pntstype == 'lght':
        gdat.numbangl = 100
        gdat.binsangl = linspace(0., gdat.maxmangl, gdat.numbangl) # [rad]
        gdat.binsanglcosi = sort(cos(gdat.binsangl))
    
    # plotting
    ## the normalized offset for text annotation of point sources in the frames
    gdat.offstext = gdat.maxmgang * 0.05
    ## figure size
    gdat.plotsize = 7
    ## size of the images
    gdat.sizeimag = 1.3 * gdat.plotsize
    ## text
    if gdat.datatype == 'mock':
        gdat.truelabl = 'Mock'
    if gdat.datatype == 'inpt':
        gdat.truelabl = gdat.strgcatl
    if gdat.strganglunit != None:
        gdat.strgxaxitotl = gdat.strgxaxi + ' [%s]' % gdat.strganglunit
        gdat.strgyaxitotl = gdat.strgyaxi + ' [%s]' % gdat.strganglunit
    else:
        gdat.strgxaxitotl = gdat.strgxaxi
        gdat.strgyaxitotl = gdat.strgyaxi

    gdat.truelablvari = gdat.truelabl + ' variable'
    gdat.truelablmiss = gdat.truelabl + ' miss'
    gdat.truelablbias = gdat.truelabl + ' off'
    gdat.truelablhits = gdat.truelabl + ' hit'
    gdat.truelablmult = gdat.truelabl + ' mult'

    gdat.truelablhost = gdat.truelabl + ' host'
    gdat.truelablsour = gdat.truelabl + ' sour'
    
    ## scaled angle axis to be plotted
    if gdat.pntstype == 'lght':
        gdat.binsanglplot = gdat.anglfact * gdat.binsangl

    # PS indices in each population
    gdat.indxpntspopl = []
    for l in gdat.indxpopl:
        gdat.indxpntspopl.append(arange(sum(gdat.maxmnumbpnts[:l]), sum(gdat.maxmnumbpnts[:l+1])))
    
    ## PSF class indices for which images will be plotted
    if gdat.numbevtt == 1:
        gdat.indxevttplot = gdat.indxevtt
    else:
        gdat.indxevttplot = concatenate((array([None]), gdat.indxevtt))
   
    if gdat.pntstype == 'lght' and gdat.pixltype != 'unbd':
        gdat.evalcirc = True
    else:
        gdat.evalcirc = False
    
    if gdat.pixltype == 'unbd':
        gdat.correxpo = False
    else:
        gdat.correxpo = True
    
    # spatial priors
    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
   
    # off-axis angle
    if gdat.modlvarioaxi or gdat.truevarioaxi:
        gdat.numboaxi = 100
        gdat.minmoaxi = 0.
        gdat.maxmoaxi = 1.1 * sqrt(2.) * gdat.maxmgangmodl
        gdat.binsoaxi = linspace(gdat.minmoaxi, gdat.maxmoaxi, gdat.numboaxi)
        gdat.binsoaxiopen = gdat.binsoaxi[:-1]
    else:
        gdat.binsoaxi = None
        gdat.numboaxi = 1
    gdat.indxoaxi = arange(gdat.numboaxi)

    # construct the true PSF
    if gdat.exprtype == 'ferm':
        retr_fermpsfn(gdat)
    if gdat.exprtype == 'chan':
        retr_chanpsfn(gdat)
    if gdat.exprtype == 'sdss':
        retr_sdsspsfn(gdat)
    if gdat.exprtype == 'hubb':
        retr_hubbpsfn(gdat)
    if gdat.exprtype == 'chem':
        gdat.truevarioaxi = False
        gdat.truepsfntype = 'singgaus'
        gdat.truepsfp = array([0.1 / gdat.anglfact])
  
    # set the PSF model to the true PSF model if a model is not specified by the user
    if gdat.modlpsfntype == None:
        gdat.modlpsfntype = gdat.truepsfntype
    if gdat.modlvarioaxi == None:
        gdat.modlvarioaxi = gdat.truevarioaxi

    # construct the PSF structure
    ## true
    gdat.numbtruepsfpform, gdat.numbtruepsfpoaxi = retr_numbpsfp(gdat.truepsfntype, gdat.truevarioaxi)
    ## model
    gdat.numbpsfpform, gdat.numbpsfpoaxi = retr_numbpsfp(gdat.modlpsfntype, gdat.modlvarioaxi)

    gdat.indxpsfpoaxi = arange(gdat.numbpsfpoaxi) 
    gdat.indxpsfpform = arange(gdat.numbpsfpform)

    gdat.numbtruepsfptotl = gdat.numbtruepsfpform + gdat.numbtruepsfpoaxi
    
    gdat.numbpsfptotl = gdat.numbpsfpform + gdat.numbpsfpoaxi
    gdat.indxpsfptotl = arange(gdat.numbpsfptotl)
    gdat.numbpsfpevtt = gdat.numbener * gdat.numbpsfptotl
   
    gdat.numbpsfptotlener = gdat.numbpsfptotl * gdat.numbener
    gdat.numbpsfp = gdat.numbpsfptotl * gdat.numbener * gdat.numbevtt
    gdat.indxpsfp = arange(gdat.numbpsfp)
    gdat.indxpsfpinit = gdat.numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)
    
    if gdat.modlvarioaxi:
        gdat.indxpsfpoaxinorm = gdat.numbpsfpform + gdat.numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)
        gdat.indxpsfpoaxiindx = gdat.numbpsfpform + gdat.numbpsfptotl * arange(gdat.numbener * gdat.numbevtt) + 1

    if gdat.truevarioaxi:
        gdat.indxtruepsfpoaxinorm = gdat.numbtruepsfpform + gdat.numbtruepsfptotl * arange(gdat.numbener * gdat.numbevtt)
        gdat.indxtruepsfpoaxiindx = gdat.numbtruepsfpform + gdat.numbtruepsfptotl * arange(gdat.numbener * gdat.numbevtt) + 1

    gdat.numbenerevtt = gdat.numbener * gdat.numbevtt

    gdat.numbbacp = gdat.numbback * gdat.numbener

    gdat.numbsigc = gdat.numbener * gdat.numbevtt





    # pixelization
    if gdat.datatype == 'mock':
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2

    gdat.numbchrototl = 5
    if gdat.pntstype == 'lght':
        gdat.numbchrollik = 7
    if gdat.pntstype == 'lens':
        gdat.numbchrollik = 3

    gdat.oaxipivt = gdat.maxmgang

    # temp
    gdat.boolintpanglcosi = False

    # number of bins
    gdat.numbbins = 10

    gdat.minmnumbpnts = zeros(gdat.numbpopl, dtype=int)
    
    # the function to measure time
    if gdat.strgfunctime == 'clck':
        gdat.functime = time.clock
    if gdat.strgfunctime == 'time':
        gdat.functime = time.time

    gdat.indxback = arange(gdat.numbback)
    
    gdat.maxmnumbpntstotl = sum(gdat.maxmnumbpnts)
    gdat.indxpntstotl = arange(gdat.maxmnumbpntstotl)
    gdat.maxmnumbpntscumr = cumsum(gdat.maxmnumbpnts)
    gdat.maxmnumbpntscuml = concatenate((array([0]), gdat.maxmnumbpntscumr[:-1]))
    
    # axes
    ## longitude
    gdat.numblgalpntsprob = gdat.numbsidepntsprob
    gdat.numbbgalpntsprob = gdat.numbsidepntsprob
    gdat.binslgalpntsprob = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbsidepntsprob + 1)
    gdat.binsbgalpntsprob = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbsidepntsprob + 1)
    gdat.indxlgalpntsprob = arange(gdat.numblgalpntsprob)
    gdat.indxbgalpntsprob = arange(gdat.numbbgalpntsprob)

    gdat.binslgal, gdat.meanlgal, gdat.difflgal, gdat.numblgal, gdat.indxlgal = tdpy.util.retr_axis(gdat.minmlgal, gdat.maxmlgal, 10)
    gdat.binsbgal, gdat.meanbgal, gdat.diffbgal, gdat.numbbgal, gdat.indxbgal = tdpy.util.retr_axis(gdat.minmbgal, gdat.maxmbgal, 10)

    # convenience variables
    gdat.indxfluxprox = arange(gdat.numbfluxprox)
    gdat.binsfluxprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
    gdat.meanfluxprox = sqrt(gdat.binsfluxprox[1:] * gdat.binsfluxprox[:-1])
    
    # spectral model
    ## total number of spectral parameters allowed
    gdat.numbspeptotl = 3
    gdat.indxspeptotl = arange(gdat.numbspeptotl)
    ## number of model spectral parameters for each population
    gdat.numbspep = retr_numbspep(gdat.spectype)
    gdat.indxspep = [arange(gdat.numbspep[l]) for l in gdat.indxpopl]
    ## plotting
    ### number of bins for histogram plots of spectral parameters
    gdat.numbspepbins = 20
    ### number of standard deviations away from the mean of spectral parameters to plot
    gdat.numbstdvspepdist = 2.
    ### minima and maxima for spectral parameters
    if gdat.datatype == 'mock':
        gdat.minmsind = min(amin(gdat.sinddistmean), amin(gdat.mocksinddistmean)) - gdat.numbstdvspepdist * max(amax(gdat.sinddiststdv), amax(gdat.mocksinddiststdv))
        gdat.maxmsind = max(amax(gdat.sinddistmean), amax(gdat.mocksinddistmean)) + gdat.numbstdvspepdist * max(amax(gdat.sinddiststdv), amax(gdat.mocksinddiststdv))
        gdat.minmcurv = min(amin(gdat.curvdistmean), amin(gdat.mockcurvdistmean)) - gdat.numbstdvspepdist * max(amax(gdat.curvdiststdv), amax(gdat.mockcurvdiststdv))
        gdat.maxmcurv = max(amax(gdat.curvdistmean), amax(gdat.mockcurvdistmean)) + gdat.numbstdvspepdist * max(amax(gdat.curvdiststdv), amax(gdat.mockcurvdiststdv))
    else:
        gdat.minmsind = gdat.sinddistmean - gdat.numbstdvspepdist * gdat.sinddiststdv
        gdat.maxmsind = gdat.sinddistmean + gdat.numbstdvspepdist * gdat.sinddiststdv
        gdat.minmcurv = gdat.curvdistmean - gdat.numbstdvspepdist * gdat.curvdiststdv
        gdat.maxmcurv = gdat.curvdistmean + gdat.numbstdvspepdist * gdat.curvdiststdv
        
    gdat.minmspep = empty(gdat.numbspeptotl)
    gdat.minmspep[0] = gdat.minmsind
    gdat.minmspep[1] = gdat.minmcurv
    gdat.minmspep[2] = gdat.minmflux
    gdat.maxmspep = empty(gdat.numbspeptotl)
    gdat.maxmspep[0] = gdat.maxmsind
    gdat.maxmspep[1] = gdat.maxmcurv
    gdat.maxmspep[2] = gdat.maxmflux
    
    gdat.binsspep = empty((gdat.numbspepbins + 1, gdat.numbspeptotl))
    for p in gdat.indxspeptotl:
        gdat.binsspep[:, p] = linspace(gdat.minmspep[p], gdat.maxmspep[p], gdat.numbspepbins + 1)
    
    ## radial
    gdat.numbgang = 10
    gdat.binsgang = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbgang + 1)

    ## azimuthal
    gdat.numbaang = 10
    gdat.binsaang = linspace(0., 2. * pi, gdat.numbaang + 1)

    # input data
    if gdat.datatype == 'inpt':
        
        path = gdat.pathdata + 'inpt/' + gdat.strgexpr
        gdat.exprdataflux = pf.getdata(path)
        
        if gdat.pixltype == 'heal':
            if gdat.exprdataflux.ndim != 3:
                raise Exception('exprdataflux should be a 3D numpy array if pixelization is HealPix.')
        else:
            if gdat.exprdataflux.ndim != 4:
                raise Exception('exprdataflux should be a 4D numpy array if pixelization is Cartesian.')
        
        if gdat.pixltype == 'cart':
            gdat.numbsidecart = gdat.exprdataflux.shape[1]
            gdat.exprdataflux = gdat.exprdataflux.reshape((gdat.exprdataflux.shape[0], gdat.numbsidecart**2, gdat.exprdataflux.shape[3]))

        gdat.numbenerfull = gdat.exprdataflux.shape[0]
        gdat.numbpixlfull = gdat.exprdataflux.shape[1]
        gdat.numbevttfull = gdat.exprdataflux.shape[2]
        gdat.indxenerfull = arange(gdat.numbenerfull)
        gdat.indxevttfull = arange(gdat.numbevttfull)
        
        if gdat.pixltype == 'heal':
            gdat.numbsideheal = int(sqrt(gdat.numbpixlfull / 12))
    
    if gdat.pixltype == 'unbd':
        gdat.indxdatasamp = arange(gdat.numbdatasamp)
        gdat.indxpixlfull = arange(gdat.numbdatasamp)
        gdat.indxpixlrofi = arange(gdat.numbdatasamp)
        gdat.apix = (2. * gdat.maxmgang)**2
    else:
        if gdat.pixltype == 'cart':
            gdat.binslgalcart = linspace(gdat.minmlgal, gdat.maxmlgal, gdat.numbsidecart + 1)
            gdat.binsbgalcart = linspace(gdat.minmbgal, gdat.maxmbgal, gdat.numbsidecart + 1)
            gdat.lgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
            gdat.bgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
            gdat.apix = (2. * gdat.maxmgang / gdat.numbsidecart)**2
            gdat.indxpixlrofi = arange(gdat.numbsidecart**2)
            indxsidecart = arange(gdat.numbsidecart)
            temp = meshgrid(indxsidecart, indxsidecart, indexing='ij')
            gdat.bgalgrid = gdat.bgalcart[temp[1].flatten()]
            gdat.lgalgrid = gdat.lgalcart[temp[0].flatten()]
        if gdat.pixltype == 'heal':
            lgalheal, bgalheal, gdat.numbpixlfull, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)
            lgalheal = deg2rad(lgalheal)
            bgalheal = deg2rad(bgalheal)
   
            gdat.indxpixlrofi = where((fabs(lgalheal) < gdat.maxmgang) & (fabs(bgalheal) < gdat.maxmgang))[0]
            
            gdat.indxpixlrofimargextd = where((fabs(lgalheal) < 1.2 * gdat.maxmgangmodl) & (fabs(bgalheal) < 1.2 * gdat.maxmgangmodl))[0]
            gdat.indxpixlrofimarg = where((fabs(lgalheal) < gdat.maxmgangmodl) & (fabs(bgalheal) < gdat.maxmgangmodl))[0]

            gdat.lgalgrid = lgalheal
            gdat.bgalgrid = bgalheal

        gdat.indxpixlfull = arange(gdat.numbpixlfull)

    # minimum angular distance from the center of the ROI
    gdat.minmgang = 1e-3
    
    if gdat.evttbins:
        # PSF class string
        gdat.strgevtt = []
        for m in gdat.indxevtt:
            gdat.strgevtt.append('PSF%d' % gdat.indxevttincl[m])
    
    # number of components
    gdat.numbcomp = 3 + zeros(gdat.numbpopl, dtype=int)
    gdat.numbcompcolr = 3 + zeros(gdat.numbpopl, dtype=int)
    if gdat.numbener > 1:
        gdat.numbcomp += gdat.numbspep + gdat.numbener - 1
        gdat.numbcompcolr += gdat.numbspep
    gdat.indxcomp = []
    for l in gdat.indxpopl:
        gdat.indxcomp.append(arange(gdat.numbcomp[l]))

    gdat.jcbnsplt = 2.**(2 - gdat.numbener)
   
    
    # component indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcomplbhl = arange(2)
    gdat.indxcompspec = 2 + gdat.indxener
    gdat.indxcompflux = 2 + gdat.indxenerfluxdist
    gdat.indxcompsind = 2 + gdat.numbener
    gdat.indxcompcurv = 2 + gdat.numbener + 1
    gdat.indxcompexpo = 2 + gdat.numbener + 1
    gdat.indxcompspep = []
    gdat.indxcompcolr = []
    gdat.indxcompunsd = []
    gdat.indxauxipara = []
    gdat.indxcomp = []
    for l in gdat.indxpopl:
        if gdat.numbener > 1:
            gdat.indxcompspep.append(arange(2 + gdat.numbener, 2 + gdat.numbener + gdat.numbspep[l]))
            gdat.indxcompcolr.append(concatenate((gdat.indxcomplbhl, gdat.indxcompflux, gdat.indxcompspep[l])))
        else:   
            gdat.indxcompcolr.append(concatenate((gdat.indxcomplbhl, gdat.indxcompflux)))
        gdat.indxcomp.append(arange(gdat.numbcomp[l]))
        gdat.indxcompunsd.append(setdiff1d(gdat.indxcomp[l], gdat.indxcompcolr[l]))
        gdat.indxauxipara.append(arange(gdat.numbcompcolr[l]))
    gdat.indxpnts = []
    for l in gdat.indxpopl:
        gdat.indxpnts.append(arange(gdat.maxmnumbpnts[l]))

    # convenience factors for CDF and ICDF transforms
    ## mean number of PS
    gdat.factmeanpnts = log(gdat.maxmmeanpnts / gdat.minmmeanpnts)
    
    ## background parameters
    gdat.factbacp = log(gdat.maxmbacp / gdat.minmbacp)
    
    ## PS parameters
    gdat.factgang = log(gdat.maxmgangmodl / gdat.minmgang)
    gdat.factflux = log(gdat.maxmflux / gdat.minmflux)
    if gdat.enerbins:
        gdat.factener = log(gdat.maxmener / gdat.minmener)

    gdat.factfluxdistslop = arctan(gdat.maxmfluxdistslop) - arctan(gdat.minmfluxdistslop)
    gdat.factfluxdistbrek = log(gdat.maxmfluxdistbrek / gdat.minmfluxdistbrek)
    gdat.factfluxdistsloplowr = arctan(gdat.maxmfluxdistsloplowr) - arctan(gdat.minmfluxdistsloplowr)
    gdat.factfluxdistslopuppr = arctan(gdat.maxmfluxdistslopuppr) - arctan(gdat.minmfluxdistslopuppr)
    
    # exposure
    if gdat.correxpo:
        if isinstance(gdat.strgexpo, float):
            if gdat.datatype == 'mock':
                if gdat.pixltype == 'heal':
                    gdat.expo = gdat.strgexpo * ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                if gdat.pixltype == 'cart':
                    gdat.expo = gdat.strgexpo * ones((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
            if gdat.datatype == 'inpt':
                gdat.expo = gdat.strgexpo * ones_like(gdat.exprdataflux)
        else:
            path = gdat.pathdata + 'inpt/' + gdat.strgexpo
            gdat.expo = pf.getdata(path)
            if amin(gdat.expo) == amax(gdat.expo):
                raise Exception('Bad input exposure map.')
                return
            if gdat.pixltype == 'cart':
                gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))

    # backgrounds
    gdat.backflux = []
    for c in gdat.indxback:
        if isinstance(gdat.strgback[c], float):
            if gdat.datatype == 'mock':
                if gdat.pixltype == 'heal':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + gdat.strgback[c]
                if gdat.pixltype == 'cart':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + gdat.strgback[c]
                if gdat.pixltype == 'unbd':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbdatasamp, gdat.numbevttfull)) + gdat.strgback[c]
            if gdat.datatype == 'inpt':
                backfluxtemp = zeros_like(gdat.exprdataflux) + gdat.strgback[c]
        else:
            path = gdat.pathdata + 'inpt/' + gdat.strgback[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'cart':
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        gdat.backflux.append(backfluxtemp)
    
    for c in gdat.indxback:
        if amin(gdat.backflux[c]) <= 0.:
            raise Exception('Background templates must be positive.')

    # only include desired energy and PSF class bins 
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    if gdat.datatype == 'inpt':
        gdat.exprdataflux = gdat.exprdataflux[gdat.indxcubeincl]
    if gdat.correxpo:
        gdat.expo = gdat.expo[gdat.indxcubeincl]
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcubeincl]
    
    # exclude voxels with vanishing exposure
    if gdat.correxpo:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.indxpixlrofi = intersect1d(gdat.indxpixlrofi, where(gdat.expo[i, :, m] > 0.)[0])
    gdat.indxcuberofi = meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.indxcube = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
     
    if gdat.pixltype != 'unbd':
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
            if gdat.verbtype > 0 and boolinitsetp:
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
        gdat.exprdataflux = gdat.exprdataflux[gdat.indxcuberofi]
   
    if gdat.correxpo:
        gdat.expofull = copy(gdat.expo)
        gdat.expo = gdat.expo[gdat.indxcuberofi]
    
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcuberofi]

    if gdat.correxpo:
        gdat.backcnts = []
        gdat.backcntstotl = zeros_like(gdat.expo)
        for c in gdat.indxback:
            backcntstemp = gdat.backflux[c] * gdat.expo * gdat.diffener[:, None, None] * gdat.apix
            gdat.backcnts.append(backcntstemp)
            gdat.backcntstotl[:] += backcntstemp 

    if gdat.pntstype == 'lght':
        if gdat.truevarioaxi:
            gdat.truefactoaxi = retr_factoaxi(gdat, gdat.binsoaxi, gdat.truepsfp[gdat.indxtruepsfpoaxinorm], gdat.truepsfp[gdat.indxtruepsfpoaxiindx])
        gdat.truepsfn = retr_psfn(gdat, gdat.truepsfp, gdat.indxener, gdat.binsangl, gdat.truepsfntype, gdat.binsoaxi, gdat.truevarioaxi)
        gdat.truefwhm = 2. * retr_psfnwdth(gdat, gdat.truepsfn, 0.5)
        
        gdat.limsangl = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        gdat.limspsfn = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.truevarioaxi:
                    psfn = gdat.truepsfn[i, :, m, 0]
                else:
                    psfn = gdat.truepsfn[i, :, m]
                maxmpsfn = amax(psfn)
                gdat.limsangl[i][m] = [0., gdat.binsangl[amax(where(psfn > 1e-6 * maxmpsfn)[0])] * gdat.anglfact]
                gdat.limspsfn[i][m] = [maxmpsfn * 1e-6, maxmpsfn]
            
    if gdat.pixltype == 'unbd':
        gdat.truepsfncdfn = roll(cumsum(gdat.truepsfn, axis=1), 1)[0, :, 0]
        gdat.truepsfncdfn[0] = 0.
        gdat.truepsfncdfn /= amax(gdat.truepsfncdfn)
        gdat.truepsfnicdfintp = interp1d(gdat.truepsfncdfn, gdat.binsangl)

    if gdat.numbpixl * gdat.maxmnumbpntstotl < 1e6 and gdat.pixltype != 'unbd':
        gdat.calcerrr = True
    else:
        gdat.calcerrr = False

    if gdat.evalcirc:
        # determine the maximum angle at which the PS flux map will be computed
        gdat.maxmangleval = empty(gdat.numbfluxprox)
        for h in gdat.indxfluxprox:
            if gdat.specfraceval == 0:
                gdat.maxmangleval[h] = 3. * gdat.maxmgangmodl
            else:   
                frac = gdat.specfraceval * gdat.binsfluxprox[0] / gdat.binsfluxprox[h+1]
                psfnwdth = retr_psfnwdth(gdat, gdat.truepsfn, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathdata + 'indxpixlprox/'
        os.system('mkdir -p %s' % path)
        path += 'indxpixlprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), 1e2 * amax(gdat.maxmangleval), gdat.numbfluxprox)
        if gdat.verbtype > 0 and boolinitsetp:
            print 'PSF evaluation will be performed up to %.3g %s for the largest flux.' % (amax(gdat.maxmangleval) * gdat.anglfact, gdat.strganglunittext)
        if os.path.isfile(path):
            if gdat.verbtype > 0 and boolinitsetp:
                print 'Previously computed nearby pixel look-up table will be used.'
                print 'Reading %s...' % path
            fobj = open(path, 'rb')
            gdat.indxpixlprox = cPickle.load(fobj)
            fobj.close()
        else:
            if gdat.verbtype > 0 and boolinitsetp:
                print 'Computing the look-up table...'
            gdat.indxpixlprox = [[] for h in range(gdat.numbfluxprox)]
            cntrsave = -1.
            # temp
            for j in gdat.indxpixl:
                if gdat.indxpixl.size == 1500**2:
                    for h in range(gdat.numbfluxprox):
                        gdat.indxpixlprox[h].append(array([0]))
                else: 
                    dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
                    dist[j] = 0.
                    for h in range(gdat.numbfluxprox):
                        indxpixlproxtemp = where(dist < gdat.maxmangleval[h])[0]
                        gdat.indxpixlprox[h].append(indxpixlproxtemp)
                cntrsave = tdpy.util.show_prog(j, gdat.indxpixl.size, cntrsave)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()


def setpfinl(gdat, boolinitsetp=False):

    gdat.maxmnumbcomp = gdat.maxmnumbpnts * gdat.numbcomp
    gdat.maxmnumbcompcumr = cumsum(gdat.maxmnumbcomp)
    gdat.maxmnumbcompcuml = concatenate((array([0]), gdat.maxmnumbcompcumr[:-1]))
    gdat.maxmnumbcomptotl = sum(gdat.maxmnumbcomp)
    
    # number of samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc

    # run tag
    gdat.rtag = retr_rtag(gdat)
    
    # plot paths
    if gdat.makeplot:
        gdat.pathplot = gdat.pathimag + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
        gdat.pathfram = gdat.pathplot + 'fram/'
        gdat.pathpost = gdat.pathplot + 'post/'
        gdat.pathinit = gdat.pathplot + 'init/'
        os.system('mkdir -p %s %s %s %s' % (gdat.pathplot, gdat.pathfram, gdat.pathpost, gdat.pathinit))
        gdat.pathdiag = gdat.pathplot + 'diag/'
        os.system('mkdir -p %s' % gdat.pathdiag)
        if gdat.optiprop:
            gdat.pathopti = gdat.pathplot + 'opti/'
            os.system('mkdir -p %s' % gdat.pathopti)
 
    # get the experimental catalog
    if gdat.exprinfo:
        
        gdat.exprcnts = None
        if gdat.exprtype == 'ferm':
            retr_fermdata(gdat)
        if gdat.exprtype == 'chan':
            retr_chandata(gdat)
    
        # rotate PS coordinates to the ROI center
        if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
            rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), gdat.anglcatlrttr], deg=True, eulertype='ZYX')
            gdat.exprbgal, gdat.exprlgal = rttr(pi / 2. - gdat.exprbgal, gdat.exprlgal)
            gdat.exprbgal = pi / 2. - gdat.exprbgal
     
        # select PSs in the ROI
        gdat.indxpntsrofi = arange(gdat.exprlgal.size, dtype=int)
        # temp
        #gdat.indxpntsrofi = intersect1d(where((fabs(gdat.exprlgal) < gdat.maxmgangmodl) & \
        #                                                    (fabs(gdat.exprbgal) < gdat.maxmgangmodl) & (amin(gdat.exprspec[0, :, :], 0) > 0.))[0], gdat.indxpntsrofi)
        
        gdat.indxpntsrofi = intersect1d(where((fabs(gdat.exprlgal) < gdat.maxmgangmodl) & \
                                                            (fabs(gdat.exprbgal) < gdat.maxmgangmodl))[0], gdat.indxpntsrofi)

        gdat.exprlgal = gdat.exprlgal[gdat.indxpntsrofi]
        gdat.exprbgal = gdat.exprbgal[gdat.indxpntsrofi]
        gdat.exprspec = gdat.exprspec[:, :, gdat.indxpntsrofi]
        if gdat.exprspep != None:
            gdat.exprspep = gdat.exprspep[gdat.indxpntsrofi, :]
        if gdat.exprcnts != None:
            gdat.exprcnts = gdat.exprcnts[:, gdat.indxpntsrofi, :]
        gdat.exprnumbpnts = gdat.exprlgal.size
    
        # compute the catalog counts based on the exposure
        gdat.exprcntscalc = empty((gdat.numbener, gdat.exprnumbpnts, gdat.numbevtt))
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                indxpixltemp = retr_indxpixl(gdat, gdat.exprbgal, gdat.exprlgal)
                gdat.exprcntscalc[i, :, m] = gdat.exprspec[0, i, :] * gdat.expo[i, indxpixltemp, m] * gdat.diffener[i]
       
        if gdat.exprcnts != None and gdat.exprlgal.size > 0 and gdat.verbtype > 0:
            if amax(fabs((gdat.exprcnts - gdat.exprcntscalc) / gdat.exprcnts)) > 0.01:
                print 'Experimental information on PS counts is inconsistent.'
        
        gdat.exprcnts = gdat.exprcntscalc

        gdat.exprgang = retr_gang(gdat.exprlgal, gdat.exprbgal)
        gdat.expraang = retr_aang(gdat.exprlgal, gdat.exprbgal)
        
        if not isfinite(gdat.exprspec).all():
            raise Exception('exprspec is not finite.')
        
        if gdat.exprnumbpnts > 0:
            gdat.exprfluxbrgt, gdat.exprfluxbrgtassc = retr_fluxbrgt(gdat, gdat.exprlgal, gdat.exprbgal, gdat.exprspec[0, gdat.indxenerfluxdist[0], :])

    # get count data
    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = gdat.exprdataflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]
    
    # load mock catalog into the reference catalog data structure
    if gdat.trueinfo:
        if gdat.datatype == 'mock':
            gdat.truemeanpnts = gdat.mocknumbpnts
        
            gdat.truelgal = []
            gdat.truebgal = []
            for l in gdat.indxpopl:
                gdat.truelgal.append(gdat.mocklgal[l])
                gdat.truebgal.append(gdat.mockbgal[l])
            
            if gdat.numbener > 1:
                gdat.truespep = []
                for l in gdat.indxpopl:
                    gdat.truespep.append(gdat.mockspep[l])
                    
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
          
            gdat.truebacp = gdat.mockbacp
            gdat.datacnts = gdat.mockdatacnts
            gdat.truefluxdistslop = gdat.mockfluxdistslop
            gdat.truefluxdistbrek = gdat.mockfluxdistbrek
            gdat.truefluxdistsloplowr = gdat.mockfluxdistsloplowr
            gdat.truefluxdistslopuppr = gdat.mockfluxdistslopuppr
    
 
    if gdat.pixltype == 'cart':
        # temp
        gdat.indxxaximaxm, gdat.indxyaximaxm = tdpy.util.retr_indximagmaxm(gdat.datacnts[0, :, 0].reshape((gdat.numbsidecart, gdat.numbsidecart)))

    if gdat.pixltype == 'unbd':
        gdat.bgalgrid = gdat.datacnts[0, :, 0, 0]
        gdat.lgalgrid = gdat.datacnts[0, :, 0, 1]
    
    # spatially averaged background flux 
    gdat.backfluxmean = zeros((gdat.numbback, gdat.numbener))
    for c in gdat.indxback:
        for i in gdat.indxener:
            if gdat.correxpo:
                gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, :] * gdat.expo[i, :, :]) / sum(gdat.expo[i, :, :])
            else:
                gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, 0])
    
    # factors in the prior expression
    gdat.priofactmeanpnts = log(1. / (log(gdat.maxmmeanpnts) - log(gdat.minmmeanpnts)))
    gdat.priofactlgalbgal = 2. * log(1. / 2. / gdat.maxmgang)
    gdat.priofactfluxdistslop = gdat.numbener * log(1. / (arctan(gdat.maxmfluxdistslop) - arctan(gdat.minmfluxdistslop)))
    gdat.priofactsplt = -2. * log(2. * gdat.maxmgangmodl) + log(gdat.radispmr) + log(2. * pi)
    # temp -- brok terms are missing

    # determine proposal probabilities
    gdat.minmlgalmarg = -gdat.maxmgangmodl
    gdat.maxmlgalmarg = gdat.maxmgangmodl
    gdat.minmbgalmarg = -gdat.maxmgangmodl
    gdat.maxmbgalmarg = gdat.maxmgangmodl

    # construct lists of possible changes to the number of PS for each PS model and the associated probabilities
    gdat.listnumbpntsmodi = []
    gdat.probnumbpntsmodi = []
    for k in gdat.indxpntstotl:
        gdat.listnumbpntsmodi.append(arange(1, k + 1))
        gdat.probnumbpntsmodi.append(1. / gdat.listnumbpntsmodi[k])
        gdat.probnumbpntsmodi[k] /= sum(gdat.probnumbpntsmodi[k])
   
    if gdat.verbtype > 1 and boolinitsetp:
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
    gdat.alphpnts = 0.4
    gdat.alphmaps = 1.
    
    ## upper limit of histograms
    gdat.limshist = [0.5, 10**ceil(log10(gdat.maxmnumbpntstotl))]

    ## marker size
    gdat.minmmrkrsize = 50
    gdat.maxmmrkrsize = 500
    ## ROI
    gdat.exttrofi = array([gdat.minmlgal, gdat.maxmlgal, gdat.minmbgal, gdat.maxmbgal])
    gdat.exttrofi *= gdat.anglfact 
    gdat.frambndrdata = gdat.maxmgang * gdat.anglfact
    gdat.frambndrmodl = gdat.maxmgangmodl * gdat.anglfact
    
    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(1000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)

    if gdat.correxpo:
        # limits on counts, which are used to bin or overplot PS counts 
        gdat.minmcnts = gdat.minmflux * mean(mean(gdat.expo, 1), 1) * gdat.diffener * gdat.factspecener
        gdat.maxmcnts = gdat.maxmflux * mean(mean(gdat.expo, 1), 1) * gdat.diffener * gdat.factspecener
        gdat.binscnts = zeros((gdat.numbener, gdat.numbfluxplot + 1))
        for i in gdat.indxener:
            gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbfluxplot + 1) # [1]
       
    ## Real data
    # true data
    if gdat.trueinfo:
        if gdat.datatype == 'inpt':
            gdat.truenumbpnts = None
            gdat.truemeanpnts = None
            gdat.truebacp = None
            
            gdat.truenumbpnts = array([gdat.exprnumbpnts], dtype=int)
            gdat.truelgal = [gdat.exprlgal]
            gdat.truebgal = [gdat.exprbgal]
            gdat.truespec = [gdat.exprspec]
            gdat.truecnts = [gdat.exprcnts]
            gdat.truespep = [gdat.exprspep]
            #gdat.truestrg = [gdat.exprstrg]
            #gdat.truestrgclss = [gdat.exprstrgclss]
            #gdat.truestrgassc = [gdat.exprstrgassc]

            gdat.trueminmflux = amin(gdat.truespec[0][0, gdat.indxenerfluxdist[0], :])
            gdat.truemaxmflux = amax(gdat.truespec[0][0, gdat.indxenerfluxdist[0], :])
            for l in gdat.indxpopl: 
                gdat.trueminmflux = min(gdat.trueminmflux, amin(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]))
                gdat.truemaxmflux = max(gdat.truemaxmflux, amax(gdat.truespec[l][0, gdat.indxenerfluxdist[0], :]))
            
            gdat.truefluxdistslop = None
            gdat.truefluxdistbrek = None
            gdat.truefluxdistsloplowr = None
            gdat.truefluxdistslopuppr = None
    
    if gdat.trueinfo and gdat.correxpo and gdat.pntstype == 'lght':
        truebackcnts = []
        gdat.truesigm = []
        for l in gdat.indxpopl:
            indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
            truebackcntstemp = zeros((gdat.numbener, gdat.truenumbpnts[l], gdat.numbevtt))
            for k in range(gdat.truebgal[l].size):
                if gdat.truevarioaxi:
                    indxoaxitemp = retr_indxoaxipnts(gdat, gdat.truelgal[l][k], gdat.truebgal[l][k])
                    fwhmtemp = gdat.truefwhm[:, :, indxoaxitemp]
                else:
                    fwhmtemp = gdat.truefwhm
                for c in gdat.indxback:
                    truebackcntstemp[:, k, :] += gdat.backflux[c][:, indxpixltemp[k], :] * gdat.expo[:, indxpixltemp[k], :] * gdat.diffener[:, None] * pi * fwhmtemp**2 / 4.
            truebackcnts.append(truebackcntstemp)
            gdat.truesigm.append(gdat.truecnts[l] / sqrt(truebackcntstemp))
    
        for l in gdat.indxpopl:
            if not isfinite(gdat.truespec[l]).all():
                print 'truespec'
                print gdat.truespec
                raise Exception('True PS parameters are not finite.')

    gdat.truefluxbrgt, gdat.truefluxbrgtassc = retr_fluxbrgt(gdat, concatenate(gdat.truelgal), concatenate(gdat.truebgal), \
                                                                                                        concatenate(gdat.truespec)[0, gdat.indxenerfluxdist[0], :])
    
    if gdat.trueinfo:
        if gdat.datatype == 'mock':
            gdat.trueminmflux = gdat.mockminmflux
            gdat.truemaxmflux = gdat.mockmaxmflux

    # set sample vector indices
    cntr = tdpy.util.cntr()
    gdat.indxfixpnumbpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxfixpmeanpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxfixpfluxdistslop = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxfixpfluxdistbrek = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxfixpfluxdistsloplowr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxfixpfluxdistslopuppr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxfixpfluxdist = concatenate((gdat.indxfixpfluxdistslop, gdat.indxfixpfluxdistbrek, gdat.indxfixpfluxdistsloplowr, gdat.indxfixpfluxdistslopuppr)) 
    gdat.indxfixphypr = concatenate((gdat.indxfixpmeanpnts, gdat.indxfixpfluxdist))

    gdat.indxfixpsigc = []
    gdat.indxfixpsigt = []
    gdat.indxfixpgamc = []
    gdat.indxfixpgamt = []
    gdat.indxfixppsff = []
    gdat.indxfixpoaxinorm = []
    gdat.indxfixpoaxiindx = []
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if gdat.modlpsfntype == 'singgaus':
                gdat.indxfixpsigc.append(cntr.incr())
            if gdat.modlpsfntype == 'singking':
                gdat.indxfixpsigc.append(cntr.incr())
                gdat.indxfixpgamc.append(cntr.incr())
            if gdat.modlpsfntype == 'doubgaus':
                gdat.indxfixpsigc.append(cntr.incr())
                gdat.indxfixpsigt.append(cntr.incr())
                gdat.indxfixppsff.append(cntr.incr())
            if gdat.modlpsfntype == 'gausking':
                gdat.indxfixpsigc.append(cntr.incr())
                gdat.indxfixpgamc.append(cntr.incr())
                gdat.indxfixpsigt.append(cntr.incr())
                gdat.indxfixppsff.append(cntr.incr())
            if gdat.modlpsfntype == 'doubking':
                gdat.indxfixpsigc.append(cntr.incr())
                gdat.indxfixpgamc.append(cntr.incr())
                gdat.indxfixpsigt.append(cntr.incr())
                gdat.indxfixpgamt.append(cntr.incr())
                gdat.indxfixppsff.append(cntr.incr())
            if gdat.modlpsfntype == 'neww':
                pass
            if gdat.modlvarioaxi:
                gdat.indxfixpoaxinorm.append(cntr.incr())
                gdat.indxfixpoaxiindx.append(cntr.incr())
    
    gdat.indxfixpsigm = gdat.indxfixpsigc + gdat.indxfixpsigt
    gdat.indxfixpgamm = gdat.indxfixpgamc + gdat.indxfixpgamt
    gdat.indxfixppsfp = gdat.indxfixpsigc + gdat.indxfixpsigt + gdat.indxfixpgamc + gdat.indxfixpgamt + gdat.indxfixppsff + gdat.indxfixpoaxinorm + gdat.indxfixpoaxiindx
    
    gdat.indxfixpbacp = []
    if gdat.pntstype == 'lght':
        gdat.indxfixpbacp = arange(gdat.numbbacp).reshape((gdat.numbback, gdat.numbener)) + cntr.incr(gdat.numbbacp)

    gdat.indxfixplenp = []
    gdat.indxfixpanglsour = []
    gdat.indxfixpanglhost = []
    gdat.indxfixpangllens = []
    if gdat.pntstype == 'lens':
        gdat.numbfixplens = 13
        gdat.indxfixplgalsour = cntr.incr()
        gdat.indxfixpbgalsour = cntr.incr()
        gdat.indxfixpfluxsour = cntr.incr()
        gdat.indxfixpsizesour = cntr.incr()
        gdat.indxfixpratisour = cntr.incr()
        gdat.indxfixpanglsour = cntr.incr()
        gdat.indxfixplgalhost = cntr.incr()
        gdat.indxfixpbgalhost = cntr.incr()
        gdat.indxfixpellphost = cntr.incr()
        gdat.indxfixpanglhost = cntr.incr()
        gdat.indxfixpsherhost = cntr.incr()
        gdat.indxfixpsanghost = cntr.incr()
        gdat.indxfixpbeinhost = cntr.incr()
        gdat.indxfixpsour = [gdat.indxfixplgalsour, gdat.indxfixpbgalsour, gdat.indxfixpfluxsour, gdat.indxfixpsizesour, gdat.indxfixpratisour, gdat.indxfixpanglsour]
        gdat.indxfixphost = [gdat.indxfixplgalhost, gdat.indxfixpbgalhost, gdat.indxfixpellphost, gdat.indxfixpanglhost, gdat.indxfixpsherhost, \
                                                                                                                                gdat.indxfixpsanghost, gdat.indxfixpbeinhost]
        gdat.indxfixplenp = gdat.indxfixpsour + gdat.indxfixphost
    
    # number of fixed-dimension parameters
    gdat.numbfixp = cntr.incr(0)
    gdat.numbfixpprim = cntr.incr(0) - gdat.numbpopl
    # indices of fixed-dimension parameters
    gdat.indxfixp = arange(gdat.numbfixp)

    # construct the fixed parameter structure
    gdat.scalfixp = zeros(gdat.numbfixp, dtype=object)
    gdat.minmfixp = zeros(gdat.numbfixp)
    gdat.maxmfixp = zeros(gdat.numbfixp)
    gdat.strgfixp = zeros(gdat.numbfixp, dtype=object)
    gdat.strgfixpunit = zeros(gdat.numbfixp, dtype=object)
    gdat.namefixp = zeros(gdat.numbfixp, dtype=object)
    gdat.factfixp = zeros(gdat.numbfixp)
    gdat.meanfixp = zeros(gdat.numbfixp)
    gdat.stdvfixp = zeros(gdat.numbfixp)
    gdat.truefixp = [None for k in gdat.indxfixp]
    gdat.cdfnminmfixp = empty(gdat.numbfixp)
    gdat.cdfndifffixp = empty(gdat.numbfixp)

    for k in gdat.indxfixp:
        if k in gdat.indxfixpnumbpnts or k in gdat.indxfixphypr:
            
            l = gdat.indxfixpnumbpnts % k
            
            if k in gdat.indxfixpnumbpnts:
                gdat.scalfixp[k] = 'pois'
                gdat.minmfixp[k] = 0
                gdat.maxmfixp[k] = gdat.maxmnumbpnts[l]
                gdat.strgfixp[k] = '$N$'
                gdat.truefixp[k] = gdat.truenumbpnts[l]
                gdat.namefixp[k] = 'numbpntspop%d' % l
                
            if k in gdat.indxfixpmeanpnts:
                gdat.scalfixp[k] = 'logt'
                gdat.minmfixp[k] = gdat.minmmeanpnts[l]
                gdat.maxmfixp[k] = gdat.maxmmeanpnts[l]
                gdat.strgfixp[k] = r'$\mu$'
                gdat.truefixp[k] = gdat.truemeanpnts[l]
                gdat.namefixp[k] = 'meanpntspop%d' % l
    
            if k == gdat.indxfixpfluxdistslop:
                gdat.scalfixp[k] = 'atan'
                gdat.minmfixp[k] = gdat.minmfluxdistslop[l]
                gdat.maxmfixp[k] = gdat.maxmfluxdistslop[l]
                gdat.strgfixp[k] = r'$\alpha$'
                gdat.truefixp[k] = gdat.truefluxdistslop[l]
                gdat.namefixp[k] = 'fluxdistsloppop%d' % l
       
            if k == gdat.indxfixpfluxdistbrek:
                gdat.scalfixp[k] = 'logt'
                gdat.minmfixp[k] = gdat.minmflux
                gdat.maxmfixp[k] = gdat.maxmflux
                gdat.strgfixp[k] = '$f_b$'
                gdat.truefixp[k] = gdat.truefluxdistbrek[l]
                gdat.namefixp[k] = 'fluxdistbrekpop%d' % l
    
            if k == gdat.indxfixpfluxdistsloplowr:
                gdat.scalfixp[k] = 'atan'
                gdat.minmfixp[k] = gdat.minmfluxdistsloplowr[l]
                gdat.maxmfixp[k] = gdat.maxmfluxdistsloplowr[l]
                gdat.strgfixp[k] = r'$\alpha_l$'
                gdat.truefixp[k] = gdat.truefluxdistsloplowr[l]
                gdat.namefixp[k] = 'fluxdistsloplowrpop%d' % l
    
            if k == gdat.indxfixpfluxdistslopuppr:
                gdat.scalfixp[k] = 'atan'
                gdat.minmfixp[k] = gdat.minmfluxdistslopuppr[l]
                gdat.maxmfixp[k] = gdat.maxmfluxdistslopuppr[l]
                gdat.strgfixp[k] = r'$\alpha_u$'
                gdat.truefixp[k] = gdat.truefluxdistslopuppr[l]
                gdat.namefixp[k] = 'fluxdistslopupprpop%d' % l

        if k in gdat.indxfixppsfp:
            if gdat.psfninfoprio:
                gdat.scalfixp[k] = 'gaus'
                n = k - gdat.indxfixppsfp[0]
                gdat.meanfixp[k] = gdat.truepsfp[n]
                gdat.stdvfixp[k] = gdat.truepsfp[n] * 0.1
            else:
                if k in gdat.indxfixpsigm:
                    gdat.scalfixp[k] = 'logt'
                    gdat.minmfixp[k] = gdat.minmsigm
                    gdat.maxmfixp[k] = gdat.maxmsigm
                if k in gdat.indxfixpgamm:
                    gdat.scalfixp[k] = 'atan'
                    gdat.minmfixp[k] = gdat.minmgamm
                    gdat.maxmfixp[k] = gdat.maxmgamm
                if k in gdat.indxfixppsff:
                    gdat.scalfixp[k] = 'atan'
                    gdat.minmfixp[k] = gdat.minmpsff
                    gdat.maxmfixp[k] = gdat.maxmpsff
                if k in gdat.indxfixpoaxinorm:
                    gdat.scalfixp[k] = 'logt'
                    gdat.minmfixp[k] = gdat.minmoaxinorm
                    gdat.maxmfixp[k] = gdat.maxmoaxinorm
                if k in gdat.indxfixpoaxiindx:
                    gdat.scalfixp[k] = 'atan'
                    gdat.minmfixp[k] = gdat.minmoaxiindx
                    gdat.maxmfixp[k] = gdat.maxmoaxiindx
                
            # strings for PSF parameters
            if k in gdat.indxfixpsigm:
                strgvarb = '\sigma'
                strgname = 'sigm'
            if k in gdat.indxfixpgamm:
                strgvarb = '\gamma'
                strgname = 'gamm'
            if k in gdat.indxfixppsff:
                strgvarb = 'f'
                strgname = 'psff'
            if k in gdat.indxfixpoaxinorm:
                strgvarb = 'a'
                strgname = 'oaxinorm'
            if k in gdat.indxfixpoaxiindx:
                strgvarb = 'b'
                strgname = 'oaxiindx'
            if gdat.modlpsfntype == 'doubgaus' or gdat.modlpsfntype == 'doubking' or gdat.modlpsfntype == 'gausking':
                if k in gdat.indxfixpsigc or k in gdat.indxfixpgamc:
                    strgcomp = 'c'
                    strgname += 'core'
                else:
                    strgcomp = 't'
                    strgname += 'tail'
            else:
                strgcomp = ''
            if gdat.numbener > 1:
                indxenertemp = gdat.indxenerincl[((k - gdat.indxfixppsfp[0]) % (gdat.numbener * gdat.numbpsfptotl)) // gdat.numbpsfptotl]
                strgenertemp = '%s' % indxenertemp
            else:
                strgenertemp = ''
            if gdat.numbevtt > 1:
                indxevtttemp = gdat.indxevttincl[(k - gdat.indxfixppsfp[0]) // (gdat.numbener * gdat.numbpsfptotl)]
                strgevtttemp = '%s' % indxevtttemp
            else:
                strgevtttemp = ''
            gdat.strgfixp[k] = r'$%s^{%s}_{%s%s}$' % (strgvarb, strgcomp, strgenertemp, strgevtttemp)
            gdat.namefixp[k] = '%spop%d_%s%s' % (strgname, l, strgenertemp, strgevtttemp)
    
        if k in gdat.indxfixpbacp:
            c = (k - gdat.indxfixpbacp[0, 0]) % gdat.numbback
            if gdat.numbener > 1:
                i = (k - gdat.indxfixpbacp[0, 0]) // gdat.numbback
                strgener = '%d' % i
            else:
                strgener = ''

            if gdat.numbback > 1:
                strgback = '%d' % c
            else:
                strgback = ''
            gdat.scalfixp[k] = 'logt'
            gdat.minmfixp[k] = gdat.minmbacp[c]
            gdat.maxmfixp[k] = gdat.maxmbacp[c]
            gdat.strgfixp[k] = '$A_{%s%s}$' % (strgback, strgener)
            gdat.namefixp[k] = 'lgalsour'
        
        if gdat.pntstype == 'lens':
            if k in gdat.indxfixplenp:
                if k == gdat.indxfixplgalsour:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = gdat.minmlgal
                    gdat.maxmfixp[k] = gdat.maxmlgal
                    gdat.strgfixp[k] = '$l_{src}$'
                    gdat.namefixp[k] = 'lgalsour'
                if k == gdat.indxfixpbgalsour:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = gdat.minmbgal
                    gdat.maxmfixp[k] = gdat.maxmbgal
                    gdat.strgfixp[k] = '$b_{src}$'
                    gdat.namefixp[k] = 'bgalsour'
                if k == gdat.indxfixpfluxsour:
                    gdat.scalfixp[k] = 'logt'
                    gdat.minmfixp[k] = gdat.minmfluxsour
                    gdat.maxmfixp[k] = gdat.maxmfluxsour
                    gdat.strgfixp[k] = '$f_{src}$'
                    gdat.namefixp[k] = 'fluxsour'
                if k == gdat.indxfixpsizesour:
                    gdat.scalfixp[k] = 'logt'
                    gdat.minmfixp[k] = gdat.minmsizesour
                    gdat.maxmfixp[k] = gdat.maxmsizesour
                    gdat.strgfixp[k] = '$a_{src}$'
                    gdat.namefixp[k] = 'sizesour'
                if k == gdat.indxfixpratisour:
                    gdat.scalfixp[k] = 'logt'
                    gdat.minmfixp[k] = gdat.minmratisour
                    gdat.maxmfixp[k] = gdat.maxmratisour
                    gdat.strgfixp[k] = '$r_{src}$'
                    gdat.namefixp[k] = 'ratisour'
                if k == gdat.indxfixpanglsour:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = 0.
                    gdat.maxmfixp[k] = 180.
                    gdat.strgfixp[k] = r'$\phi_{src}$'
                    gdat.namefixp[k] = 'anglsour'
                if k == gdat.indxfixplgalhost:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = gdat.minmlgal
                    gdat.maxmfixp[k] = gdat.maxmlgal
                    gdat.strgfixp[k] = '$l_{host}$'
                    gdat.namefixp[k] = 'lgalhost'
                if k == gdat.indxfixpbgalhost:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = gdat.minmbgal
                    gdat.maxmfixp[k] = gdat.maxmbgal
                    gdat.strgfixp[k] = '$b_{host}$'
                    gdat.namefixp[k] = 'bgalhost'
                if k == gdat.indxfixpellphost:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = 0.
                    gdat.maxmfixp[k] = 0.5
                    gdat.strgfixp[k] = r'$\epsilon_{host}$'
                    gdat.namefixp[k] = 'ellphost'
                if k == gdat.indxfixpanglhost:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = 0.
                    gdat.maxmfixp[k] = 180.
                    gdat.strgfixp[k] = r'$\phi_{host}$'
                    gdat.namefixp[k] = 'anglhost'
                if k == gdat.indxfixpsherhost:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = gdat.minmsherhost
                    gdat.maxmfixp[k] = gdat.maxmsherhost
                    gdat.strgfixp[k] = r'$\gamma_e$'
                    gdat.namefixp[k] = 'sherhost'
                if k == gdat.indxfixpsanghost:
                    gdat.scalfixp[k] = 'self'
                    gdat.minmfixp[k] = 0.
                    gdat.maxmfixp[k] = 180.
                    gdat.strgfixp[k] = r'$\phi_{\gamma}$'
                    gdat.namefixp[k] = 'sanghost'
                if k == gdat.indxfixpbeinhost:
                    gdat.scalfixp[k] = 'logt'
                    gdat.meanfixp[k] = gdat.minmbeinhost
                    gdat.stdvfixp[k] = gdat.maxmbeinhost
                    gdat.strgfixp[k] = '$R_{host}$'
                    gdat.namefixp[k] = 'beinhost'
        
        if gdat.scalfixp[k] == 'self':
            gdat.factfixp[k] = gdat.maxmfixp[k] - gdat.minmfixp[k]
        if gdat.scalfixp[k] == 'logt':
            gdat.factfixp[k] = log(gdat.maxmfixp[k] / gdat.minmfixp[k])
        if gdat.scalfixp[k] == 'atan':
            gdat.factfixp[k] = arctan(gdat.maxmfixp[k]) - arctan(gdat.minmfixp[k])
        if gdat.scalfixp[k] == 'gaus':
            gdat.minmfixp[k] = gdat.meanfixp[k] - 3. * gdat.stdvfixp[k]
            gdat.maxmfixp[k] = gdat.meanfixp[k] + 3. * gdat.stdvfixp[k]
        if gdat.scalfixp[k] == 'eerr':
            gdat.cdfnminmfixp[k], gdat.cdfndifffixp[k] = retr_eerrnorm(gdat.minmfixp[k], gdat.maxmfixp[k], gdat.meanfixp[k], gdat.stdvfixp[k])
        
        if k in gdat.indxfixpfluxdistbrek:
            gdat.strgfixpunit[k] = gdat.strgfixp[k] + ' [%s]' % gdat.strgfluxunit
        elif k in gdat.indxfixpsigc or k in gdat.indxfixpsigt:
            gdat.strgfixpunit[k] = gdat.strgfixp[k] + ' [%s]' % gdat.strganglunit
        else:
            gdat.strgfixpunit[k] = gdat.strgfixp[k]

    gdat.numbpara = gdat.numbfixp + gdat.maxmnumbcomptotl
    gdat.indxsampcomp = arange(gdat.indxfixp[-1] + 1, gdat.numbpara)
    gdat.indxpara = arange(gdat.numbpara)
    
    if gdat.verbtype > 1 and boolinitsetp:
        print 'gdat.numbpara'
        print gdat.numbpara
        print 'gdat.numbfixp'
        print gdat.numbfixp
        print 'gdat.maxmnumbcomp'
        print gdat.maxmnumbcomp
        print 'gdat.maxmnumbpnts'
        print gdat.maxmnumbpnts
        print 'gdat.maxmnumbcomptotl'
        print gdat.maxmnumbcomptotl
        print 'gdat.maxmnumbpntstotl'
        print gdat.maxmnumbpntstotl
        print 'gdat.strgfixp'
        print gdat.strgfixp
        print 'gdat.scalfixp'
        print gdat.scalfixp
        print

    # transdimensional parameters
    gdat.numbtrap = gdat.numbpara - gdat.indxsampcomp[0]
    gdat.indxsamptrap = arange(gdat.indxsampcomp[0], gdat.numbpara)
   
    # proposals
    gdat.strgprop = []
    cntr = tdpy.util.cntr()
    gdat.indxfixpused = []
    
    gdat.indxfixpused.append(gdat.indxfixpnumbpnts)
    
    ## mean number of PS
    gdat.indxpropmeanpnts = cntr.incr()
    gdat.strgprop.append('meanpnts')
    gdat.indxfixpused.append(gdat.indxfixpmeanpnts)

    ## flux distribution
    gdat.indxpropfluxdist = cntr.incr()
    gdat.strgprop.append('fluxdistslop')
    for l in gdat.indxpopl:
        if gdat.fluxdisttype[l] == 'powr':  
            gdat.indxfixpused.append(array([gdat.indxfixpfluxdistslop[l]]))
        if gdat.fluxdisttype[l] == 'brok':  
            gdat.indxfixpused.append(array([gdat.indxfixpfluxdistbrek[l]]))
            gdat.indxfixpused.append(array([gdat.indxfixpfluxdistsloplowr[l]]))
            gdat.indxfixpused.append(array([gdat.indxfixpfluxdistslopuppr[l]]))

    # PSF parameters
    gdat.indxproppsfp = None
    if gdat.proppsfp:
        gdat.indxproppsfp = cntr.incr()
        gdat.strgprop.append('psfp')
        gdat.indxfixpused.append(gdat.indxfixppsfp)
    
    gdat.indxpropbacp = None
    # background parameters
    if gdat.propbacp:
        gdat.indxpropbacp = cntr.incr()
        gdat.strgprop.append('bacp')
        gdat.indxfixpused.append(gdat.indxfixpbacp.flatten())
    
    gdat.indxproplenp = None
    # lensing parameters
    if gdat.pntstype == 'lens':
        gdat.indxproplenp = cntr.incr()
        gdat.strgprop.append('lens')
        gdat.indxfixpused.append(gdat.indxfixplenp)
    
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
    
    # spep
    gdat.indxpropspep = None
    if gdat.numbener > 1:
        gdat.strgprop.append('spep')
        gdat.indxpropspep = cntr.incr()

    gdat.numbprop = len(gdat.strgprop)
    gdat.indxprop = arange(gdat.numbprop)
            
    gdat.indxproplpri = [gdat.indxpropmeanpnts, gdat.indxpropfluxdist]
    gdat.indxpropllik = setdiff1d(setdiff1d(gdat.indxprop, gdat.indxpropmeanpnts), gdat.indxpropfluxdist)

    gdat.probproptran = 0.4
    gdat.probpropbrth = 1.
    gdat.stdvproppara = 0.1 + zeros(gdat.numbfixp)
    gdat.indxpropbrde = array([gdat.indxpropbrth, gdat.indxpropdeth])
    gdat.indxpropspme = array([gdat.indxpropsplt, gdat.indxpropmerg])
    gdat.probtranmaxm = array([0., 1.])
    gdat.probtranminm = array([1., 0.])
    gdat.probtran = array([.5, .5])

    gdat.indxfixpused = concatenate(gdat.indxfixpused)
    gdat.indxfixpusedprim = gdat.indxfixpused[gdat.numbpopl:]
    

    
    if gdat.modlvarioaxi:
        gdat.indxfixppsfpoaxinorm = gdat.indxfixppsfp[0] + gdat.indxpsfpoaxinorm
        gdat.indxfixppsfpoaxiindx = gdat.indxfixppsfp[0] + gdat.indxpsfpoaxiindx
        gdat.indxfixppsfpoaxi = sort(concatenate((gdat.indxfixppsfpoaxinorm, gdat.indxfixppsfpoaxiindx)))

    listindxsampunsd = []
    numbsampcumu = 0
    for l in gdat.indxpopl:
        for k in gdat.indxpnts[l]:
            listindxsampunsd.append(gdat.indxsampcomp[0] + numbsampcumu + k * gdat.numbcomp[l] + gdat.indxcompunsd[l])
        numbsampcumu += gdat.maxmnumbpnts[l] * gdat.numbcomp[l]
    gdat.indxsampunsd = concatenate(listindxsampunsd)
    
    ## flux
    gdat.minmfluxplot = gdat.minmflux
    gdat.maxmfluxplot = gdat.maxmflux
    if gdat.trueinfo:
        gdat.minmfluxplot = min(gdat.minmfluxplot, gdat.trueminmflux)
        gdat.maxmfluxplot = max(gdat.maxmfluxplot, gdat.truemaxmflux)
    gdat.binsfluxplot = logspace(log10(gdat.minmfluxplot), log10(gdat.maxmfluxplot), gdat.numbfluxplot + 1)
    gdat.meanfluxplot = sqrt(gdat.binsfluxplot[1:] * gdat.binsfluxplot[:-1])
    gdat.difffluxplot = gdat.binsfluxplot[1:] - gdat.binsfluxplot[:-1]
    
    ## color
    gdat.sindcdfnminm, gdat.sindcdfndiff = retr_eerrnorm(gdat.minmsind, gdat.maxmsind, gdat.sinddistmean, gdat.sinddiststdv)
    if gdat.datatype == 'mock':
        gdat.mocksindcdfnminm, gdat.mocksindcdfndiff = retr_eerrnorm(gdat.minmsind, gdat.maxmsind, gdat.mocksinddistmean, gdat.mocksinddiststdv)

    gdat.binssind = linspace(gdat.minmsind, gdat.maxmsind, gdat.numbspepbins + 1)
    gdat.meansind = (gdat.binssind[1:] + gdat.binssind[:-1]) / 2.
    gdat.diffsind = gdat.binssind[1:] - gdat.binssind[:-1]

    gdat.minmspecplot = gdat.minmfluxplot * gdat.factspecener
    gdat.maxmspecplot = gdat.maxmfluxplot * gdat.factspecener
    gdat.binsspecplot = gdat.binsfluxplot[None, :] * gdat.factspecener[:, None]
    gdat.meanspecplot = empty((gdat.numbener, gdat.numbfluxplot))
    for i in gdat.indxener:
        gdat.meanspecplot[i, :] = sqrt(gdat.binsspecplot[i, 1:] * gdat.binsspecplot[i, :-1])

    # determine the indices of true point sources, which will be compared againts the model sources
    if gdat.trueinfo:
        gdat.indxtruepntscomp = []
        for l in gdat.indxpopl:
            indxtruepntstemp = where((fabs(gdat.truelgal[l]) < gdat.maxmgangcomp) & (fabs(gdat.truebgal[l]) < gdat.maxmgangcomp))[0]
            gdat.indxtruepntscomp.append(indxtruepntstemp)

    # sanity checks
    # temp
    if (fabs(gdat.datacnts - rint(gdat.datacnts)) > 1e-3).any() and boolinitsetp:
        print 'Fractional counts!'

    if amin(gdat.datacnts) < 0. and boolinitsetp:
        print 'Negative counts!'

    # plotting
    numbtickcbar = 10
    gdat.tickdatacnts = empty((gdat.numbener, numbtickcbar))
    gdat.labldatacnts = empty((gdat.numbener, numbtickcbar), dtype=object)
    if gdat.pixltype != 'unbd':
        gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix / gdat.diffener
    else:
        gdat.datafluxmean = array([gdat.numbdatasamp / gdat.apix])
        
    if gdat.pixltype != 'unbd':
        gdat.datacntsmean = mean(sum(gdat.datacnts, 2), 1)
    
    if gdat.pixltype != 'unbd':
        gdat.minmdatacnts = amin(amin(gdat.datacnts, 1), 1)
        if gdat.satumaps:
            gdat.maxmdatacnts = ceil((amax(sum(gdat.datacnts, 2), 1) - gdat.datacntsmean) * 0.05 + gdat.datacntsmean)
        else:
            gdat.maxmdatacnts = amax(sum(gdat.datacnts, 2), 1)
    else:
        gdat.minmdatacnts = array([gdat.numbdatasamp / gdat.apix]) * 1e-1
        gdat.maxmdatacnts = array([gdat.numbdatasamp / gdat.apix]) * 1e1
    
    if gdat.pixltype != 'unbd':
        gdat.maxmresicnts = ceil(gdat.maxmdatacnts * 0.1)
        gdat.tickresicnts = empty((gdat.numbener, numbtickcbar + 1))
        gdat.lablresicnts = empty((gdat.numbener, numbtickcbar + 1), dtype=object)
        if gdat.calcerrr:
            gdat.maxmerrrcnts = ceil(gdat.maxmdatacnts * 0.02)
            gdat.maxmerrr = ones(gdat.numbener) 
            gdat.tickerrrcnts = empty((gdat.numbener, numbtickcbar + 1))
            gdat.lablerrrcnts = empty((gdat.numbener, numbtickcbar + 1), dtype=object)
            gdat.tickerrr = empty((gdat.numbener, numbtickcbar + 1))
            gdat.lablerrr = empty((gdat.numbener, numbtickcbar + 1), dtype=object)
        if gdat.scalmaps == 'asnh':
            gdat.minmdatacnts = arcsinh(gdat.minmdatacnts)
            gdat.maxmdatacnts = arcsinh(gdat.maxmdatacnts)
            gdat.maxmresicnts = arcsinh(gdat.maxmresicnts)
        for i in gdat.indxener:
            gdat.tickdatacnts[i, :] = linspace(gdat.minmdatacnts[i], gdat.maxmdatacnts[i], numbtickcbar)
            gdat.tickresicnts[i, :] = linspace(-gdat.maxmresicnts[i], gdat.maxmresicnts[i], numbtickcbar + 1)
            
            if gdat.calcerrr:
                gdat.tickerrrcnts[i, :] = linspace(-gdat.maxmerrrcnts[i], gdat.maxmerrrcnts[i], numbtickcbar + 1)
                gdat.tickerrr[i, :] = linspace(-gdat.maxmerrr[i], gdat.maxmerrr[i], numbtickcbar + 1)
            for k in range(numbtickcbar +1):
                gdat.lablresicnts[i, k] = '%.3g' % gdat.tickresicnts[i, k]
            
                if gdat.calcerrr:
                    gdat.lablerrrcnts[i, k] = '%.3g' % gdat.tickerrrcnts[i, k]
                    gdat.lablerrr[i, k] = '%.3g' % gdat.tickerrr[i, k]
                if gdat.scalmaps == 'asnh':
                    gdat.lablresicnts[i, k] = '%.3g' % sinh(gdat.tickresicnts[i, k])
                else:
                    gdat.lablresicnts[i, k] = '%.3g' % gdat.tickresicnts[i, k]
                if k != numbtickcbar:
                    if gdat.scalmaps == 'asnh':
                        gdat.labldatacnts[i, k] = '%.3g' % sinh(gdat.tickdatacnts[i, k])
                    else:
                        gdat.labldatacnts[i, k] = '%.3g' % gdat.tickdatacnts[i, k]
         
    if gdat.verbtype > 1 and boolinitsetp:
        if gdat.pntstype == 'lght' and gdat.pixltype != 'unbd':
            print 'Memory budget: indxpixlprox'
            totl = 0.
            for h in gdat.indxfluxprox:
                for n in gdat.indxpixl:
                    totl += sys.getsizeof(gdat.indxpixlprox[h][n]) / 2.**20
            print '%.4g MB' % totl


def retr_fluxbrgt(gdat, lgal, bgal, flux):

    if lgal.size == 0:
        fluxbrgt = array([0.])
        fluxbrgtassc = array([0.])
    else:
        indxbrgt = argmax(flux)
        fluxbrgt = flux[indxbrgt]
        dir1 = array([lgal, bgal])#[:, None]
        dir2 = array([lgal[indxbrgt], bgal[indxbrgt]])
        dist = retr_angldist(gdat, dir1, dir2)
        indxbrgtassc = where(dist < gdat.anglassc)[0]
        fluxbrgtassc = flux[indxbrgtassc]
        fluxbrgt = repeat(fluxbrgt, fluxbrgtassc.size)

    return fluxbrgt, fluxbrgtassc


def retr_indxoaxipnts(gdat, lgal, bgal):

    dir1 = array([lgal, bgal])[:, None]
    oaxi = retr_angldist(gdat, dir1, array([0., 0.]))
    # temp -- check that digitize works as expected
    indxoaxipnts = digitize(oaxi[0], gdat.binsoaxiopen)

    return indxoaxipnts


def init_figr(gdat, strgplot, gdatmodi=None, indxenerplot=None, indxevttplot=None, indxpoplplot=None, pathfold=None):

    if pathfold == None:
        pathfold = gdat.pathfram

    figr, axis = plt.subplots(figsize=(gdat.sizeimag, gdat.sizeimag))
    
    if indxenerplot == None:
        strgener = ''
    else:
        strgener = '%d' % gdat.indxenerincl[indxenerplot]
    
    if indxevttplot == None:
        strgevtt = 'A'
    else:
        strgevtt = '%d' % gdat.indxevttincl[indxevttplot]
    
    if indxpoplplot == None:
        strgpopl = ''
    else:
        strgpopl = '%d' % indxpoplplot

    if gdatmodi == None:
        strgswep = ''
    else:
        strgswep = '_swep%09d' % gdatmodi.cntrswep
    
    path = '%s%s%s%s%s%s.pdf' % (pathfold, strgplot, strgener, strgevtt, strgpopl, strgswep)
   
    axis.set_xlabel(gdat.strgxaxitotl)
    axis.set_ylabel(gdat.strgyaxitotl)
    if indxenerplot != None and gdat.numbener > 1 or indxevttplot != None and gdat.numbevtt > 1:
        if indxenerplot != None and gdat.numbener > 1:
            titl = gdat.strgbinsener[indxenerplot]
        else:
            titl += ', ' + gdat.strgevtt[indxevttplot]
        axis.set_title(titl)

    return figr, axis, path


def draw_frambndr(gdat, axis):
    
    outr = max(gdat.frambndrmodl, gdat.frambndrdata)
    axis.set_xlim([-outr, outr])
    axis.set_ylim([-outr, outr])
    innr = min(gdat.frambndrmodl, gdat.frambndrdata)
    axis.axvline(innr, ls='--', alpha=gdat.alphmrkr, color='black')
    axis.axvline(-innr, ls='--', alpha=gdat.alphmrkr, color='black')
    axis.axhline(innr, ls='--', alpha=gdat.alphmrkr, color='black')
    axis.axhline(-innr, ls='--', alpha=gdat.alphmrkr, color='black')


def retr_scat(gdat, axis, maps, thisindxener, thisindxevtt):

    draw_frambndr(gdat, axis)
    
    scat = axis.scatter(maps[thisindxener, :, thisindxevtt, 0], maps[thisindxener, :, thisindxevtt, 1], alpha=gdat.alphmaps, facecolor='black', s=5)

    return scat


def retr_imag(gdat, axis, maps, thisindxener, thisindxevtt, cmap='Reds', mean=False, vmin=None, vmax=None, scal=None):

    if vmin == None and vmax != None:
        vmin = -vmax
    
    draw_frambndr(gdat, axis)
   
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
    
    # rescale the map
    if scal == 'asnh':
        maps = arcsinh(maps)
    
    imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest', vmin=vmin, vmax=vmax, alpha=gdat.alphmaps)
    
    return imag


def make_cbar(gdat, axis, imag, indxenerplot=None, tick=None, labl=None):

    # make a color bar
    if indxenerplot != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)
        if tick != None and labl != None:
            cbar.set_ticks(tick)
            cbar.set_ticklabels(labl)
    return cbar


def make_catllabl(gdat, axis):

    axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                label='Sample', marker='+', linewidth=2, color='b')
    if gdat.trueinfo:
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                    label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                    label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                    label=gdat.truelablhits, marker='*', linewidth=2, color='g')
    axis.legend(bbox_to_anchor=[0.5, 1.1], loc='center', ncol=2)
        

def supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot, trueonly=False):

    # true catalog
    if gdat.trueinfo:
        
        ## get the true catalog
        mrkrsize = retr_mrkrsize(gdat, gdat.truespec[indxpoplplot][0, gdat.indxenerfluxdist, :].flatten())
        lgal = copy(gdat.truelgal[indxpoplplot])
        bgal = copy(gdat.truebgal[indxpoplplot])
        numbpnts = int(gdat.truenumbpnts[indxpoplplot])
        if not trueonly:
            
            ## associations
            ### missed
            indx = gdatmodi.indxtruepntsassc[indxpoplplot].miss
            axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablmiss, marker='*', linewidth=2, color='g')
            
            ### biased
            indx = gdatmodi.indxtruepntsassc[indxpoplplot].bias[indxenerplot]
            axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                        label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
            
            ### hit
            indx = gdatmodi.indxtruepntsassc[indxpoplplot].hits[indxenerplot]
            
        else:
            indx = arange(lgal.size)
        
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablhits, marker='x', linewidth=2, color='g')
        
        if gdat.pntstype == 'lens':
           
            ## host
            axis.scatter(gdat.anglfact * gdat.mocklgalhost, gdat.anglfact * gdat.mockbgalhost, alpha=gdat.alphpnts, label=gdat.truelablhits, \
                                                                                                                    s=300, marker='D', linewidth=2, color='g')
            axis.add_patch(plt.Circle((gdat.anglfact * gdat.mocklgalhost, gdat.anglfact * gdat.mockbgalhost), \
                                                                                        gdat.fluxfactplot * gdat.mockbeinhost, edgecolor='g', facecolor='none', lw=2))
            
            ## source
            axis.scatter(gdat.anglfact * gdat.mocklgalsour, gdat.anglfact * gdat.mockbgalsour, alpha=gdat.alphpnts, label=gdat.truelablhits, \
                                                                                                                           s=300, marker='>', linewidth=2, color='g')

            for k in range(lgal.size):
                axis.add_patch(plt.Circle((gdat.anglfact * lgal[k], gdat.anglfact * bgal[k]), \
                                gdat.fluxfactplot * gdat.truespec[indxpoplplot][0, gdat.indxenerfluxdist[0], k], edgecolor='g', facecolor='none', lw=2))

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
    if not trueonly:
        mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[indxpoplplot][gdat.indxenerfluxdist, :]])
        lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[indxpoplplot]]
        bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[indxpoplplot]]
        axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphpnts, label='Sample', marker='+', linewidth=2, color='b')
        if gdat.pntstype == 'lens':
            
            # subhalos
            for k in range(lgal.size):
                axis.add_artist(plt.Circle((gdat.anglfact * lgal[k], gdat.anglfact * bgal[k]), \
                        gdat.fluxfactplot * gdat.truespec[indxpoplplot][0, gdat.indxenerfluxdist[0], k], edgecolor='b', facecolor='none', ls='--', lw=2))

            ## host
            lgalhost = gdatmodi.thissampvarb[gdat.indxfixplgalhost]
            bgalhost = gdatmodi.thissampvarb[gdat.indxfixpbgalhost]
            beinhost = gdatmodi.thissampvarb[gdat.indxfixpbeinhost]
            axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, alpha=gdat.alphpnts, label='Model Host', s=300, marker='s', linewidth=2, color='b')
            axis.add_patch(plt.Circle((gdat.anglfact * lgalhost, gdat.anglfact * bgalhost), gdat.fluxfactplot * beinhost, edgecolor='b', facecolor='none', lw=2, ls='--'))
            
            ## source
            lgalsour = gdatmodi.thissampvarb[gdat.indxfixplgalsour]
            bgalsour = gdatmodi.thissampvarb[gdat.indxfixpbgalsour]
            axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, alpha=gdat.alphpnts, label='Model Source', s=300, marker='<', linewidth=2, color='b')


def retr_levi(listllik):
    
    minmlistllik = amin(listllik)
    levi = log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    return levi


def retr_info(listllik, levi):
    
    info = mean(listllik) - levi

    return info


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


def retr_angldist(gdat, dir1, dir2):
    
    if gdat.pixltype == 'cart':
        angldist = sqrt((dir1[0, :] - dir2[0])**2 + (dir1[1, :] - dir2[1])**2)
    else:
        angldist = angdist(dir1, dir2)

    return angldist


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
        dist = retr_angldist(gdat, dir1, dir2)
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
    minmflux = gdat.minmspecplot[indxenerthis]
    maxmflux = gdat.maxmspecplot[indxenerthis]

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
        

