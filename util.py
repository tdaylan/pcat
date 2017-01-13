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


def retr_spec(gdat, flux, sind, curv, expo, spectype):

    #if isscalar(flux):
    #    flux = array([flux])

    #if sind.ndim == 1:
    #    sind = sind[None, :]
    
    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if spectype == 'powr':
            spec = flux[None, :] * gdat.enernorm[:, None]**(-sind[None, :])
        if spectype == 'curv':
            spec = flux[None, :] * gdat.enernorm[:, None]**(-sind[None, :] - gdat.factlogtenerpivt[:, None] * curv[None, :])
        if spectype == 'expo':
            spec = flux[None, :] * gdat.enernorm[:, None]**(-sind[None, :]) * exp(gdat.enerexpofact[:, None] / expo[None, :])

    return spec


def retr_indx(gdat, indxpntsfull, spectype):
    
    indxsamplgal = []
    indxsampbgal = []
    indxsampflux = []
    indxsampspec = []
    indxsampsind = []
    indxsampcurv = [[] for l in gdat.indxpopl]
    indxsampexpo = [[] for l in gdat.indxpopl]
    indxsampcompcolr = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            indxsamplgaltemp = gdat.indxsampcomp[0] + gdat.maxmnumbcompcuml[l] + array(indxpntsfull[l], dtype=int) * gdat.numbcomp[l]
            indxsamplgal.append(indxsamplgaltemp)
            indxsampbgal.append(indxsamplgaltemp + 1)
            indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], gdat.numbener, 0) + repeat(gdat.indxener, len(indxpntsfull[l])).reshape(gdat.numbener, -1))
            indxsampflux.append(indxsampspec[l][gdat.indxenerfluxdist[0], :])
            if gdat.numbener > 1:
                indxsampsind.append(indxsamplgaltemp + 2 + gdat.numbener)
                if spectype[l] == 'curv':
                    indxsampcurv[l] = indxsampsind[l] + 1
                if spectype[l] == 'expo':
                    indxsampexpo[l] = indxsampsind[l] + 1

            indxsampcompcolr.append(repeat(indxsamplgaltemp, gdat.numbcompcolr[l]) + tile(gdat.indxcompcolr[l], len(indxpntsfull[l])))
             
    return indxsamplgal, indxsampbgal, indxsampflux, indxsampspec, indxsampsind, indxsampcurv, indxsampexpo, indxsampcompcolr


def retr_fluxhistprio(gdat, l, sampvarb):
    
    meanpnts = sampvarb[gdat.indxfixpmeanpnts[l]]
    if gdat.fluxdisttype == 'powr':
        fluxdistslop = sampvarb[gdat.indxfixpfluxdistslop[l]]  
    if gdat.fluxdisttype == 'brok':
        fluxdistbrek = sampvarb[gdat.indxsampfluxdistbrek[l]]  
        fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]  
        fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]  

    if gdat.fluxdisttype == 'powr':
        fluxhistprio = meanpnts * pdfn_flux_powr(gdat, gdat.meanfluxplot, fluxdistslop) * gdat.deltfluxplot
    if gdat.fluxdisttype == 'brok':
        fluxhistprio = meanpnts * pdfn_flux_brok(gdat, gdat.meanfluxplot, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr) * gdat.deltfluxplot
    
    return fluxhistprio
            

def retr_plotpath(gdat, strg, gdatmodi):
    
    if gdatmodi == None:
        path = gdat.pathpost + strg + '.pdf'
    else:
        path = gdat.pathfram + strg + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


def retr_pntsflux(gdat, lgal, bgal, spec, psfnintp, varioaxi, evalcirc):
  
    if gdat.verbtype > 1:
        print 'retr_pntsflux'

    numbpnts = lgal.size
    if gdat.pixltype == 'unbd':
        pntsfluxsing = zeros((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt, 2))
    else:
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
        
        if varioaxi:
            indxoaxitemp = retr_indxoaxipnts(gdat, lgal[k], bgal[k])
            psfntemp = psfnintp[indxoaxitemp](dist)
        else:
            psfntemp = psfnintp(dist)
        
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                pntsfluxsing[k, i, indxpixltemp, m] = spec[i, k] * psfntemp[i, :, m]
                
    # sum contributions from all PS
    pntsflux = sum(pntsfluxsing, 0) 
    
    return pntsflux


def retr_mapslght(gdat, bacp, pntsflux, tempindx):
    
    modlflux = pntsflux[tempindx]
    for c in gdat.indxback:
        modlflux += bacp[c*gdat.numbener+gdat.indxener, None, None] * gdat.backflux[c][tempindx]        

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


def lpdf_gaus(para, meanpara, stdvpara):
    
    lpdf = -0.5  * log(2. * pi) * stdvpara - 0.5 * (para - meanpara)**2 / stdvpara**2
    
    return lpdf


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


def icdf_fixp(gdat, strg, fixpunit, thisindxfixp):

    scalfixp = getattr(gdat, strg + 'scalfixp')[thisindxfixp]

    if scalfixp == 'self' or scalfixp == 'logt' or scalfixp == 'atan':
        
        minmfixp = getattr(gdat, strg + 'minmfixp')[thisindxfixp]
        factfixp = getattr(gdat, strg + 'factfixp')[thisindxfixp]
        
        if scalfixp == 'self':
            fixp = icdf_self(fixpunit, minmfixp, factfixp)
        elif scalfixp == 'logt':
            fixp = icdf_logt(fixpunit, minmfixp, factfixp)
        elif scalfixp == 'atan':
            fixp = icdf_atan(fixpunit, minmfixp, factfixp)
    
    elif scalfixp == 'gaus' or scalfixp == 'eerr':
        
        meanfixp = getattr(gdat, strg + 'meanfixp')[thisindxfixp]
        stdvfixp = getattr(gdat, strg + 'stdvfixp')[thisindxfixp]
        
        if scalfixp == 'eerr':
            
            cdfnminmfixp = getattr(gdat, strg + 'cdfnminmfixp')[thisindxfixp]
            cdfndifffixp = getattr(gdat, strg + 'cdfndifffixp')[thisindxfixp]
        
            fixp = icdf_eerr(fixpunit, meanfixp, stdvfixp, cdfnminmfixp, cdfndifffixp)
        else:
            fixp = icdf_gaus(fixpunit, meanfixp, stdvfixp)

    elif scalfixp == 'pois':
        fixp = fixpunit
    else:
        raise Exception('Scaling of the parameter is unrecognized.')

    return fixp


def retr_thisindxprop(gdat, gdatmodi):

    # temp
    #gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampflux, gdatmodi.thisindxsampspec, gdatmodi.thisindxsampsind, \
    #                                                                        gdatmodi.thisindxsampcompcolr = retr_indx(gdat, gdatmodi.thisindxpntsfull)

    # initialize the Boolean flag indicating the type of transdimensional proposal
    gdatmodi.propbrth = False
    gdatmodi.propdeth = False
    gdatmodi.propsplt = False
    gdatmodi.propmerg = False
   
    # index of the population in which a transdimensional proposal will be made
    gdatmodi.indxpoplmodi = choice(gdat.indxpopl)
    if rand() < gdat.probbrde:
        
        ## births and deaths
        if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]:
            gdatmodi.propdeth = True
        elif gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.minmnumbpnts[gdatmodi.indxpoplmodi]:
            gdatmodi.propbrth = True
        else:
            if rand() < 0.5:
                gdatmodi.propbrth = True
            else:
                gdatmodi.propdeth = True
        
        if gdatmodi.propbrth:
            gdatmodi.thisindxprop = gdat.indxpropbrth
        else:
            gdatmodi.thisindxprop = gdat.indxpropdeth

    else:
        
        ## splits and merges
        if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]:
            gdatmodi.propmerg = True
        elif gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.minmnumbpnts[gdatmodi.indxpoplmodi]:
            gdatmodi.propbrth = True
        else:
            if rand() < 0.5:
                gdatmodi.propsplt = True
            else:
                gdatmodi.propmerg = True
        
        if gdatmodi.propsplt:
            gdatmodi.thisindxprop = gdat.indxpropsplt
        else:
            gdatmodi.thisindxprop = gdat.indxpropmerg

    if gdat.verbtype > 1:
        print 
        print 'retr_thisindxprop()'
        print 'propbrth'
        print gdatmodi.propbrth
        print 'propdeth'
        print gdatmodi.propdeth
        print 'propsplt'
        print gdatmodi.propsplt
        print 'propmerg'
        print gdatmodi.propmerg
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


def retr_cntsmaps(gdat, fluxmaps, cart=False):

    if cart:
        cntsmaps = fluxmaps * gdat.expocart * gdat.apix
        if gdat.enerbins:
            cntsmaps *= gdat.diffener[:, None, None, None]
    else:
        
        cntsmaps = fluxmaps * gdat.expo * gdat.apix
        if gdat.enerbins:
            cntsmaps *= gdat.diffener[:, None, None]
        
    return cntsmaps


def retr_cntsbackfwhm(gdat, bacp, fwhm):

    varioaxi = len(fwhm.shape) == 3
    cntsbackfwhm = zeros_like(fwhm)
    for c in gdat.indxback:
        indxbacp = c * gdat.numbener + gdat.indxener
        if varioaxi:
            cntsback = bacp[indxbacp, None, None, None] * gdat.backflux[c][:, :, :, None] * gdat.expo[:, :, :, None] * \
                                                                                                gdat.diffener[:, None, None, None] * pi * fwhm[:, None, :, :]**2 / 4.
        else:
            cntsback = bacp[indxbacp, None, None] * gdat.backflux[c] * gdat.expo * pi * fwhm[:, None, :]**2 / 4.
            if gdat.enerbins:
                cntsback *= gdat.diffener[:, None, None]
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


def retr_probpois(data, modl):
    
    prob = data * log(modl) - modl - sp.special.gammaln(data + 1)

    return prob

    #if gdatmodi.propsplt or gdatmodi.propmerg:
    #    
    #    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
    #        lprbfrst = log(pdfn_flux_powr(gdat, gdatmodi.fluxfrst, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]))
    #        lprbseco = log(pdfn_flux_powr(gdat, gdatmodi.fluxseco, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]))
    #        lprbpare = log(pdfn_flux_powr(gdat, gdatmodi.fluxpare, gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]]))
    #    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
    #        lprbfrst += log(pdfn_flux_brok(gdat, gdatmodi.fluxfrst, \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]], \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]))
    #        
    #        lprbseco += log(pdfn_flux_brok(gdat, gdatmodi.fluxseco, \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]], \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]))

    #        lprbpare -= log(pdfn_flux_brok(gdat, gdatmodi.fluxpare, \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistbrek[gdatmodi.indxpoplmodi]], \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistsloplowr[gdatmodi.indxpoplmodi]], \
    #                           gdatmodi.thissampvarb[gdat.indxfixpfluxdistslopuppr[gdatmodi.indxpoplmodi]]))
    #    if gdatmodi.propsplt:
    #        gdatmodi.nextlpri[gdat.indxlprispme] += lprbfrst
    #        gdatmodi.nextlpri[gdat.indxlprispme] += lprbseco
    #        gdatmodi.nextlpri[gdat.indxlprispme] -= lprbpare
    #    else:
    #        gdatmodi.nextlpri[gdat.indxlprispme] += lprbpare
    #        gdatmodi.nextlpri[gdat.indxlprispme] -= lprbfrst
    #        gdatmodi.nextlpri[gdat.indxlprispme] -= lprbseco
            
    # temp
    ## split
    # P(f1)P(l1)P(b1)P(s1)P(f2)P(l2)P(b2)P(s2) / P(f0)P(l0)P(b0)P(s0)P(uf)P(ur)P(up)P(us)
    # P(f1)P(f2)P(l2)P(b2) / P(f0)P(uf)P(ur)P(up)
    # P(f1)P(f2) / P(f0)

        
def retr_sampvarb(gdat, indxpntsfull, samp, strg):
    
    if strg == 'true':
        strgtemp = strg
    else:
        strgtemp = ''
    spectype = getattr(gdat, strgtemp + 'spectype')
    indxsamplgal, indxsampbgal, indxsampflux, indxsampspec, indxsampsind, indxsampcurv, indxsampexpo, indxsampcompcolr = retr_indx(gdat, indxpntsfull, spectype) 
    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxfixpnumbpnts] = samp[gdat.indxfixpnumbpnts]
    
    for k in gdat.indxfixp:
        sampvarb[k] = icdf_fixp(gdat, '', samp[k], k)
    
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl)
            sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgangmodl, 2. * gdat.maxmgangmodl) 
            if gdat.fluxdisttype == 'powr':
                sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_powr(samp[indxsampspec[l][gdat.indxenerfluxdist, :]], gdat.minmflux, gdat.maxmflux, \
                                                                                                                                         sampvarb[gdat.indxfixpfluxdistslop[l]])
            if gdat.fluxdisttype == 'brok':
                fluxunit = samp[indxsampspec[l][gdat.indxenerfluxdist[0], :]]
                fluxdistbrek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
                fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
                fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
                sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_brok(fluxunit, gdat.minmflux, gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
            
            if gdat.numbener > 1:
                sampvarb[indxsampsind[l]] = icdf_gaus(samp[indxsampsind[l]], sampvarb[gdat.indxfixpsinddistmean[l]], sampvarb[gdat.indxfixpsinddiststdv[l]])
                if gdat.spectype[l] == 'curv':
                    sampvarb[indxsampcurv[l]] = icdf_gaus(samp[indxsampcurv[l]], gdat.curvddistmean[l], gdat.curvdiststdv[l])
                if gdat.spectype[l] == 'expo':
                    sampvarb[indxsampexpo[l]] = icdf_logt(samp[indxsampexpo[l]], gdat.minmener, gdat.factener)
            
                print 'l'
                print l
                print 'gdat.spectype'
                print gdat.spectype
                print 'indxsampexpo'
                print indxsampexpo
                print 
                sampvarb[indxsampspec[l]] = retr_spec(gdat, sampvarb[indxsampflux[l]], sampvarb[indxsampsind[l]], \
                                                                                    curv=sampvarb[indxsampcurv[l]], expo=sampvarb[indxsampexpo[l]], spectype=gdat.spectype[l])

    return sampvarb
    

def retr_mrkrsize(gdat, flux):

    mrkrsize = (log(flux) - log(gdat.minmflux)) / (log(gdat.maxmflux) - log(gdat.minmflux)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


def retr_hubbpsfn(gdat):

    gdat.exprpsfp = array([0.05]) / gdat.anglfact
    gdat.exprvarioaxi = False


def retr_sdsspsfn(gdat):
   
    gdat.exprpsfp = array([0.25 / gdat.anglfact, 1.7e6, 1.9, 0.25 / gdat.anglfact, 2.1e6, 2.])
    gdat.exprvarioaxi = False


def retr_chanpsfn(gdat):

    gdat.exprpsfp = array([0.3 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.6e-1, 2.])
    gdat.exprvarioaxi = True
   

def retr_fermpsfn(gdat):
   
    gdat.exprvarioaxi = False
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
    gdat.exprpsfp = zeros((gdat.numbener * numbpsfpform * gdat.numbevtt))
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            gdat.exprpsfp[indxfermpsfptemp] = fermform[:, m, k]

    # calculate the scale factor
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    

def updt_samp(gdat, gdatmodi):
   
    if gdat.verbtype > 1:
        print 'updt_samp()'

    # update the sample and the unit sample vectors
    gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
    gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -2] = gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1]
    
    # update the log-prior
    gdatmodi.thislpritotl = copy(gdatmodi.nextlpritotl)
    
    # update the log-likelihood
    gdatmodi.thislliktotl = copy(gdatmodi.nextlliktotl)
        
    gdatmodi.thisindxpntsfull = deepcopy(gdatmodi.nextindxpntsfull)


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

    with open(gdat.pathinpt + 'chancatl.txt', 'r') as thisfile:
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

    path = gdat.pathinpt + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
    listhdun = ap.io.fits.open(path)
    wcso = ap.wcs.WCS(listhdun[0].header)
   
    skycobjt = ap.coordinates.SkyCoord("galactic", l=lgalchan, b=bgalchan, unit='deg')
    rascchan = skycobjt.fk5.ra.degree
    declchan = skycobjt.fk5.dec.degree

    indxpixllgal = 1490
    indxpixlbgal = 1510

    # temp 0 or 1 makes a difference!
    lgalchan, bgalchan = wcso.wcs_world2pix(rascchan, declchan, 0)
    lgalchan -= gdat.numbsidecart / 2 + indxpixllgal
    bgalchan -= gdat.numbsidecart / 2 + indxpixlbgal
    lgalchan *= gdat.sizepixl
    bgalchan *= gdat.sizepixl

    gdat.exprbgal = lgalchan
    gdat.exprlgal = bgalchan
    
    gdat.exprspec = zeros((3, gdat.numbener, gdat.exprlgal.size))
    gdat.exprcnts = zeros((gdat.numbener, gdat.exprlgal.size, gdat.numbevtt))

    gdat.exprspec[0, 0, :] = fluxchansoft * 0.624e12
    gdat.exprspec[0, 1, :] = fluxchanhard * 0.624e12 / 16.
    
    # temp
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :]

    # temp
    gdat.exprspec[where(gdat.exprspec < 0.)] = 0.

    gdat.exprcnts[0, :, 0] = cntschansoft
    gdat.exprcnts[1, :, 0] = cntschanhard

    #gdat.exprsind = -log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :]) / log(gdat.enernorm)
    
    #gdat.exprstrg = lgalstrg
    #gdat.exprstrgclss = lgalchanclss
    #gdat.exprstrgassc = lgalchanassc

    #indxsort = argsort(fluxchansoft)[::-1]
    #gdat.exprlgal = gdat.exprlgal[indxsort][:150]
    #gdat.exprbgal = gdat.exprbgal[indxsort][:150]
    #gdat.exprspec = gdat.exprspec[:, :, indxsort][:150]
    #gdat.exprcnts = gdat.exprcnts[:, indxsort][:150]


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
    

def retr_rtag(gdat):
    
    rtag = '%d' % (gdat.numbswep)
        
    return rtag


def retr_gaus(gdat, gdatmodi, indxsamp, stdv, inpl=False):
    
    if gdat.probrand > 0.:
        if rand() < gdat.probrand:
            gdatmodi.drmcsamp[indxsamp, 1] = rand()
            thisbool = False
        else:
            thisbool = True
    else:
        thisbool = True
    
    if thisbool:
        if inpl:
            indx = 1
        else:
            indx = 0
        if isinstance(stdv, float):
            gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, indx] + normal(scale=stdv)
        else:
            for k in range(stdv.size):
                gdatmodi.drmcsamp[indxsamp[k], 1] = gdatmodi.drmcsamp[indxsamp[k], indx] + normal(scale=stdv[k])
            
       
def retr_deflcutf(angl, anglscal, anglcutf):
    
    fracscal = angl / anglscal
    fraccutf = anglscal / anglcutf

    deflcutf = bein / 2. / fracscal**2 * fraccutf**2 / (fraccutf**2 + 1)**2 * ((fraccutf**2 + 1. + 2. * (fracscal**2 - 1.) ) * acos(1. / fracscal) / sqrt(fracscal**2 - 1.) + \
            pi * fraccutf + (fraccutf**2 - 1.) * log(fraccutf) + sqrt(fracscal**2 + fraccutf**2) * ((fraccutf - 1. / fraccutf) * log(sqrt(fracscal**2 + fraccutf**2) \
                        / fracscal - fraccutf) - pi))
    
    return deflcutf


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
    print '%22s %14s %14s %14s %14s %14s' % ('name', 'thissamp', 'nextsamp', 'thissampvarb', 'nextsampvarb', 'diffsampvarb')
    for k in gdat.indxpara:
        if k == gdat.numbfixp:
            print
        if k < gdat.numbfixp:
            name = gdat.namefixp[k]
        else:
            name = ''
        print '%22s %14.4g %14.4g %14.4g %14.4g %14.4g' % (name, getattr(gdatmodi, 'drmcsamp')[k, 0], getattr(gdatmodi, 'drmcsamp')[k, 1], gdatmodi.thissampvarb[k], \
                                                                                    gdatmodi.nextsampvarb[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k])
    print
    return
    for strg in ['this', 'next']:
        print '%sindxpntsfull' % strg
        print getattr(gdatmodi, strg + 'indxpntsfull')
        print '%sindxsamplgal' % strg
        print getattr(gdatmodi, strg + 'indxsamplgal')
        print '%sindxsampbgal' % strg
        print getattr(gdatmodi, strg + 'indxsampbgal')
        print '%sindxsampspec' % strg
        print getattr(gdatmodi, strg + 'indxsampspec')
        print '%sindxsampspep' % strg
        print getattr(gdatmodi, strg + 'indxsampspep')
        print '%sindxsampcompcolr' % strg
        print getattr(gdatmodi, strg + 'indxsampcompcolr')
        print


def retr_prop(gdat, gdatmodi):
 
    if gdat.verbtype > 1:
        print 'retr_prop()'
    
    for k in gdat.indxactvprop:
        retr_gaus(gdat, gdatmodi, gdat.indxfixpactvprop[k], gdatmodi.stdvstdp[k])
    
    for k in gdat.indxiact:
        gdatmodi.nextsampvarb[gdat.indxsampiact[k]] = gdatmodi.thissampvarb[gdat.indxsampiact[k]]
   
    for k in gdat.indxfixpdist:
        gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, '', gdatmodi.drmcsamp[k, -1], k)

    # rescale the unit sample vector if a hyperparameter controlling the distribution of PS properties is being updated
    for l in gdat.indxpopl: 
        
        ## flux distribution
        if gdat.fluxdisttype == 'powr':
            fluxunit = cdfn_flux_powr(gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]], gdat.minmflux, gdat.maxmflux, gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslop[l]])
        if gdat.fluxdisttype == 'brok':
            fluxunit = cdfn_flux_brok(gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]], gdat.minmflux, gdat.maxmflux, gdatmodi.nextsampvarb[gdat.indxfixpfluxdistbrek[l]], \
                                            gdatmodi.nextsampvarb[gdat.indxfixpfluxdistsloplowr[l]], gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslopuppr[l]])
        
        gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :], -1] = fluxunit
    
        ## color distribution
        if gdat.numbener > 1:
            sindunit = cdfn_gaus(gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[l][:, 0]], gdatmodi.nextsampvarb[gdat.indxfixpsinddistmean[l]], \
                                                                                                                gdatmodi.nextsampvarb[gdat.indxfixpsinddiststdv[l]])
            gdatmodi.drmcsamp[gdatmodi.thisindxsampspep[l][:, 0], -1] = sindunit

    # number of unit sample vector elements to be modified
    numbcompmodi = gdat.numbcomp[gdatmodi.indxpoplmodi]
    numbcompcolrmodi = gdat.numbcompcolr[gdatmodi.indxpoplmodi]
    
    gdatmodi.nextindxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
    
    if gdatmodi.propbrth:
        
        # find an empty slot in the PS list
        for k in range(gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]):
            if not k in gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi]:
                indxpntsbrth = k
                break

        # sample indices to add the new PS
        indxsamptran = gdat.indxsampcomp[0] + gdat.maxmnumbcompcuml[gdatmodi.indxpoplmodi] + indxpntsbrth * gdat.numbcomp[gdatmodi.indxpoplmodi] + \
                                                                                                                                 gdat.indxcompcolr[gdatmodi.indxpoplmodi]
        
        # sample auxiliary variables
        gdatmodi.auxipara = rand(numbcompcolrmodi)
        
        gdatmodi.drmcsamp[indxsamptran, -1] = gdatmodi.auxipara

        if gdat.verbtype > 1:
            print 'auxipara'
            print gdatmodi.auxipara
            print 'numbcompcolrmodi'
            print numbcompcolrmodi
            print 'numbcompmodi'
            print numbcompmodi
                
    # death
    if gdatmodi.propdeth:
        
        # occupied PS index to be killed
        dethindxindxpnts = choice(arange(gdatmodi.drmcsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi], -2], dtype=int))
        
        # PS index to be killed
        indxpntsdeth = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][dethindxindxpnts]

        if gdat.verbtype > 1:
            print 'dethindxpnts: ', indxpntsdeth
            print 'dethindxindxpnts: ', dethindxindxpnts
            print
  
    ## birth
    if gdatmodi.propbrth or gdatmodi.propsplt:
        
        # change the number of PS
        gdatmodi.drmcsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi], -1] = gdatmodi.drmcsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi], -2] + 1
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].append(indxpntsbrth)

    ## death
    if gdatmodi.propdeth or gdatmodi.propmerg:
        temp = deepcopy(gdatmodi.thisindxpntsfull)[gdatmodi.indxpoplmodi]
        indxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
        indxpntsfull[gdatmodi.indxpoplmodi].remove(indxpntsdeth)
        
        # change the number of PS
        gdatmodi.drmcsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi], -1] = gdatmodi.drmcsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi], -2] - 1
        
        # remove the PS from the occupied PS list
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].remove(indxpntsdeth)
    else:
        indxpntsfull = gdatmodi.thisindxpntsfull
    
    indxsamplgal, indxsampbgal, indxsampflux, indxsampspec, indxsampspep, indxsampcompcolr = retr_indx(gdat, indxpntsfull, gdat.spectype)
    
    # PSs
    for l in gdat.indxpopl:
        stdvlbhl = gdatmodi.stdvstdp[gdat.indxstdplgal] / (gdatmodi.thissampvarb[indxsampflux[l]] / gdat.minmflux)
        for k in range(len(indxpntsfull[l])):
            retr_gaus(gdat, gdatmodi, indxsamplgal[l][k], stdvlbhl[k])
            retr_gaus(gdat, gdatmodi, indxsampbgal[l][k], stdvlbhl[k])
        retr_gaus(gdat, gdatmodi, indxsampflux[l], gdatmodi.stdvstdp[gdat.indxstdpflux], inpl=True)
        if gdat.numbener > 1:
            retr_gaus(gdat, gdatmodi, indxsampsind[l], gdatmodi.stdvstdp[gdat.indxstdpsind], inpl=True)
            if gdat.spectype[l] == 'expo':
                retr_gaus(gdat, gdatmodi, indxsampexpo[l], gdatmodi.stdvstdp[gdat.indxstdpexpo], inpl=True)
            if gdat.spectype[l] == 'curv':
                retr_gaus(gdat, gdatmodi, indxsampcurv[l], gdatmodi.stdvstdp[gdat.indxstdpcurv], inpl=True)
   
    gdatmodi.nextsampvarb = retr_sampvarb(gdat, gdatmodi.nextindxpntsfull, gdatmodi.drmcsamp[:, 1], 'next')
   
    if gdat.numbtrap > 0:
        gdatmodi.indxsampchec = concatenate((gdat.indxfixpactvprop, concatenate(indxsampcompcolr)))
        if gdatmodi.propbrth:
            gdatmodi.indxsampmodi = concatenate((gdat.indxfixpnumbpnts, gdatmodi.indxsampchec, indxsamptran))
        else:
            gdatmodi.indxsampmodi = concatenate((gdat.indxfixpnumbpnts, gdatmodi.indxsampchec))
    else:
        gdatmodi.indxsampchec = gdat.indxfixpactvprop
        gdatmodi.indxsampmodi = gdat.indxfixpactvprop
    
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbpntsmodi = 3
        
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        gdatmodi.spltindxindxpnts = choice(thisindxindxpnts)
    
        # update the sample vector
        gdatmodi.indxsampfrst = gdat.indxsampcomp[0] + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + \
                                                int(gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]) * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp[gdatmodi.indxpoplmodi]

        gdatmodi.indxsampseco = gdat.indxsampcomp[0] + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + indxpntsbrth * gdat.numbcomp[gdatmodi.indxpoplmodi]
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
            show_samp(gdat, gdatmodi)

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

    if gdatmodi.propmerg:
        
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

    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and not gdatmodi.boolreje:
        
        ## Jacobian
        jcbnfacttemp = gdatmodi.fluxpare * fabs(gdatmodi.auxipara[1] * (sin(gdatmodi.auxipara[2]) * cos(gdatmodi.auxipara[2]) + cos(gdatmodi.auxipara[2])**2))
        if gdatmodi.propsplt:
            gdatmodi.jcbnfact = jcbnfacttemp
        else:
            gdatmodi.jcbnfact = 1. / jcbnfacttemp
        
        ## combinatorial factor
        thisnumbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]
        if gdatmodi.propsplt:
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

    factoaxi = 1. + norm[:, None, None] * (bins[None, None, :] / gdat.oaxipivt)**indx[:, None, None]
    
    return factoaxi


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, binsoaxi=None, varioaxi=None, strgpara=''):

    numbpsfpform = getattr(gdat, strgpara + 'numbpsfpform')
    numbpsfptotl = getattr(gdat, strgpara + 'numbpsfptotl')
    
    indxpsfpinit = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    if varioaxi:
        indxpsfpoaxinorm = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp]
        indxpsfpoaxiindx = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp] + 1

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
   
    if psfntype == 'singgaus':
        sigc = psfp[indxpsfpinit]
        if varioaxi:
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
        else:
            sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
        
    elif psfntype == 'singking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        psfn = retr_singking(scalangl, sigc, gamc)
        if varioaxi:
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
            gamc = gamc[:, None, :, None]
        else:
            sigc = sigc[:, None, :]
            gamc = gamc[:, None, :]
        
    elif psfntype == 'doubgaus':
        frac = psfp[indxpsfpinit]
        sigc = psfp[indxpsfpinit+1]
        sigt = psfp[indxpsfpinit+2]
        if varioaxi:
            frac = frac[:, None, :, None]
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
            sigt = sigt[:, None, :, None] * factoaxi[:, None, :, :]
        else:
            frac = frac[:, None, :]
            sigc = sigc[:, None, :]
            sigt = sigt[:, None, :]
        psfn = retr_doubgaus(scalangl, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfp[indxpsfpinit]
        sigc = psfp[indxpsfpinit+1]
        sigt = psfp[indxpsfpinit+2]
        gamt = psfp[indxpsfpinit+3]
        if varioaxi:
            frac = frac[:, None, :, None]
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
            sigt = sigt[:, None, :, None] * factoaxi[:, None, :, :]
            gamt = gamt[:, None, :, None]
        else:
            frac = frac[:, None, :]
            sigc = sigc[:, None, :]
            sigt = sigt[:, None, :]
            gamt = gamt[:, None, :]
        psfn = retr_gausking(scalangl, frac, sigc, sigt, gamt)
        
    elif psfntype == 'doubking':
        frac = psfp[indxpsfpinit]
        sigc = psfp[indxpsfpinit+1]
        gamc = psfp[indxpsfpinit+2]
        sigt = psfp[indxpsfpinit+3]
        gamt = psfp[indxpsfpinit+4]
        if varioaxi:
            frac = frac[:, None, :, None]
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
            gamc = gamc[:, None, :, None]
            sigt = sigt[:, None, :, None] * factoaxi[:, None, :, :]
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


def retr_axis(gdat, strg, minm=None, maxm=None, numb=None, bins=None, scal='self'):
    
    if bins == None:
        if scal == 'self':
            bins = linspace(minm, maxm, numb + 1)
            mean = (bins[1:] + bins[:-1]) / 2.
        else:
            bins = logspace(log10(minm), log10(maxm), numb + 1)
            mean = sqrt(bins[1:] * bins[:-1])
    else:
        if scal == 'self':
            mean = (bins[1:] + bins[:-1]) / 2.
        else:
            mean = sqrt(bins[1:] * bins[:-1])
        numb = mean.size
    indx = arange(numb)
    delt = diff(bins) 

    setattr(gdat, 'bins' + strg, bins)
    setattr(gdat, 'mean' + strg, mean)
    setattr(gdat, 'delt' + strg, delt)
    setattr(gdat, 'numb' + strg, numb)
    setattr(gdat, 'indx' + strg, indx)


def retr_unit(lgal, bgal):

    xaxi = cos(bgal) * cos(lgal)
    yaxi = -cos(bgal) * sin(lgal)
    zaxi = sin(bgal)

    return xaxi, yaxi, zaxi


def retr_psec(gdat, conv):

    # temp
    psec = (abs(fft.fft2(conv))**2)[:gdat.numbsidecart/2, :gdat.numbsidecart/2] * 1e-3

    return psec
   

def retr_psecodim(gdat, psec):

    psecodim = zeros(gdat.numbwvecodim)
    for k in gdat.indxwvecodim:
        indxwvec = where((gdat.meanwvec > gdat.binswvecodim[k]) & (gdat.meanwvec < gdat.binswvecodim[k+1]))
        psecodim[k] = mean(psec[indxwvec])
        #indxwvec = where((gdat.meanwvec > gdat.binswvecodim[k]) & (gdat.meanwvec < gdat.binswvecodim[k+1]))[0]
        #if indxwvec.size > 0:
        #    psecodim[k] = mean(psec[indxwvec])
    
    return psecodim


def retr_randunitpsfp(gdat):

    while True:
        randunitpsfp = rand(gdat.numbpsfp)
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
                    indx = m * gdat.numbpsfpevtt + i * gdat.numbpsfpform
                    thisbool = thisbool and randunitpsfp[indx+indxpar1] > randunitpsfp[indx+indxpar0]
            if thisbool:
                break

    return randunitpsfp


def retr_varb(gdat, strg, gdatmodi=None, perc='medi'):
        
    if gdatmodi != None:
        varb = gdatmodi.thissampvarb[getattr(gdat, 'indx' + strg)]
    else:
        varb = getattr(gdat, perc + strg)

    return varb


def retr_massfrombein(bein):
    
    mass = 10**12 * (bein / (1.8 * pi / 3600. / 180.))**2
    
    return mass


def retr_eerrnorm(minmvarb, maxmvarb, meanvarb, stdvvarb):
   
    cdfnminm = 0.5 * (sp.special.erf((minmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfnmaxm = 0.5 * (sp.special.erf((maxmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfndiff = cdfnmaxm - cdfnminm
    
    return cdfnminm, cdfndiff


def retr_numbspep(spectype):
    
    numbpopl = len(spectype)
    numbspep = empty(numbpopl, dtype=int)
    liststrgspep = [[] for l in range(numbpopl)]
    liststrgfluxspep = [[] for l in range(numbpopl)]
    for l in range(numbpopl):
        if spectype[l] == 'powr':
            liststrgspep[l] = ['sind']
        if spectype[l] == 'expo':
            liststrgspep[l] = ['sind', 'expo']
        if spectype[l] == 'curv':
            liststrgspep[l] = ['sind', 'curv']
        liststrgfluxspep[l] = ['flux'] + liststrgspep[l]
        numbspep[l] = len(liststrgspep[l]) 

    return numbspep, liststrgspep, liststrgfluxspep
    

def retr_detrcatl(gdat):
    
    # find the total number of elements across all the samples
    numbelem = 0
    for n in gdat.indxsamptotl:
        for l in gdat.indxpopl:
            numbelem += len(gdat.listindxpntsfull[n][l])

    indxpoplelem = empty(numbelem, dtype=int)
    indxsamptotlelem = empty(numbelem, dtype=int)

    indxpntsfull = deepcopy(gdat.listindxpntsfull)

    while len(indxpntsfull):
        
        numbassc = 0
        cros = []
        for n in gdat.indxsamptotl:
            for l in gdat.indxpopl:
                for k in range(len(gdat.listindxpntsfull[n][l])):
                    for nn in gdat.indxsamptotl:
                        for ll in gdat.indxpopl:
                            for kk in range(len(gdat.listindxpntsfull[nn][ll])):
                                
                                if l != ll or n == nn:
                                    continue
                    
                                ## compue distance

                                ## preselect elements to compute the cross correlation againts

                                ## compute derivative
                                deritemp = 0.
                                for e in gdat.indxcomp[l]:
                                    deritemp += 0.

                                if deritemp > derithrs:
                                    numbassc += 1
                                cros.append([n, l, k, nn, ll, kk, deritemp])

        # cluster elements based on their correlation
        #listindxclus =  
        listindxpnts = [arange(gdat.listlgal[l][n].size) for n in gdat.indxsamptotl]

        listlgalstck = concatenate(gdat.listlgal[l])
        listbgalstck = concatenate(gdat.listbgal[l])
        listfluxstck = concatenate(gdat.listflux[l])
        if gdat.numbener > 1:
            listspepstck = concatenate(gdat.listspep[l])
            
        indxsampstck = concatenate(listindxpnts)
        print 'indxsampstck'
        print summgene(indxsampstck)

        numbpntsstck = listlgalstck.size
        for k in range(numbpntsstck):
            indxsamptotlself = indxsampstck[k]
            print 'indxsamptotlself'
            print indxsamptotlself
            print 'gdat.listlgal[l]'
            print len(gdat.listlgal[l])
            print
            listlgalstcktemp = setdiff1d(listlgalstck, gdat.listlgal[l][indxsamptotlself])
            listbgalstcktemp = setdiff1d(listbgalstck, gdat.listbgal[l][indxsamptotlself])
            listfluxstcktemp = setdiff1d(listfluxstck, gdat.listflux[l][indxsamptotlself])

            dist = ((listlgalstck[k] - listlgalstck) / listlgalstck)**2 + ((listbgalstck[k] - listbgalstck) / listbgalstck)**2 + \
                                                                                                        ((listfluxstck[k] - listfluxstck) / listfluxstck)**2 
            #if gdat.numbener > 1:
            #    dist += ((listspepstck[k] - listspepstck) / listspepstck)**2
            
            indxpntsstck = where(dist < 0.1)[0]
            
        indxpntsstcksort = argsort(numbpntsstck)
        
        detrcatl = []
        for a in range(numbpntsstck):
            
            if len(indxpntsstck) > 0:
                
                # delete the member PS from the stacked catalog
                delete(listcatlstck, indxpntsstcksort[a])
                
                # calculate the posterior median of the PS in the group
                medicatlstck = median(listcatlstck)
                
                # append the median PS to the list
                detrcatl.append(catl) 


def retr_conv(gdat, defl):
    
    conv = (gradient(defl[:, :, 0], gdat.sizepixl, axis=0) + gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) / 2.

    return conv


def diag_gdatmodi(gdatmodi, gdatmodiprev):

    listvalu = []
    listattr = []
    print 'diag_gdatmodi'
    for attr, valu in gdatmodi.__dict__.iteritems():
       
        boolmodi = False
        
        try:    
            valuprev = getattr(gdatmodiprev, attr)
        except:
            pass
            continue

        if attr == 'thissampvarb' or attr == 'drmcsamp':
            
            indx = where(valu - valuprev != 0.)
            if indx[0].size > 0:
                print valu[indx]
            
        if isinstance(valuprev, ndarray):
            try:
                indx = where(valu - valuprev != 0.)
                if indx[0].size > 0:
                    boolmodi = True
            except:
                print 'ndarry failed'
                print 'attr'
                print attr
                print
                #boolmodi = True
                
        elif isinstance(valuprev, list):
            continue
            for k, item in enumerate(valu):
                if isinstance(item, list):
                    for l, itemitem in enumerate(item):
                        if valuprev[k] != item:
                            boolmodi = True
                    print item
                if valuprev[k] != item:
                    boolmodi = True
        elif isinstance(valuprev, (int, bool, float)):
            
            try:
                if valu != valuprev:
                    boolmodi = True
            except:
                boolmodi = True
        else:
            print 'other than numpy'
            print type(valuprev)
            print

        if boolmodi:
            #print 'attr'
            print attr
            if False:
                if isinstance(valuprev, ndarray):
                    if len(indx) > 1:
                        print 'valuprev'
                        print valuprev[indx[0]]
                        print 'valu'
                        print valu[indx[0]]
                    else:
                        print 'valuprev'
                        print valuprev
                        print 'valu'
                        print valu
                else:
                    print 'valuprev'
                    print valuprev
                    print 'valu'
                    print valu
            #print
    print


def retr_pntscnts(gdat, lgal, bgal, spec):
    
    indxpixltemp = retr_indxpixl(gdat, bgal, lgal)
    cnts = zeros((gdat.numbener, lgal.size, gdat.numbevtt))
    for k in range(lgal.size):
        cnts[:, k, :] += spec[:, k, None] * gdat.expo[:, indxpixltemp[k], :] * gdat.diffener[:, None]
    
    return cnts


def setpinit(gdat, boolinitsetp=False):

    # number of samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc

    # run tag
    gdat.rtag = retr_rtag(gdat)
    
    # paths
    ## data
    gdat.pathpixlcnvt = gdat.pathdata + 'pixlcnvt/'
    gdat.pathpixlprox = gdat.pathdata + 'pixlprox/'
    ## plot
    if gdat.makeplot:
        gdat.pathplot = gdat.pathimag + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
        gdat.pathinit = gdat.pathplot + 'init/'
        gdat.pathdiag = gdat.pathplot + 'diag/'
        gdat.pathfram = gdat.pathplot + 'fram/'
        gdat.pathpost = gdat.pathplot + 'post/'
        gdat.pathpostfixp = gdat.pathpost + 'fixp/'
        for strg in ['lpri', 'llik']:
            setattr(gdat, 'pathpostdelt%s' % strg, gdat.pathpost + 'delt%s/' % strg)
            setattr(gdat, 'pathpostdelt%saccp' % strg, gdat.pathpost + 'delt%saccp/' % strg)
            setattr(gdat, 'pathpostdelt%sreje' % strg, gdat.pathpost + 'delt%sreje/' % strg)
        if gdat.probtran > 0. and gdat.probbrde < 1.:
            gdat.pathpostspmr = gdat.pathpost + 'spmr/'
        if gdat.optiprop:
            gdat.pathopti = gdat.pathplot + 'opti/'
        if gdat.makeanim:
            gdat.pathanim = gdat.pathplot + 'anim/'
    
    for attr, valu in gdat.__dict__.iteritems():
        if attr.startswith('path'):
            os.system('mkdir -p %s' % valu)
 
    # spectral model
    ## total number of spectral parameters allowed
    gdat.numbspeptotl = 3
    gdat.indxspeptotl = arange(gdat.numbspeptotl)
    ## number of model spectral parameters for each population
    gdat.numbspep, gdat.liststrgspep, gdat.liststrgfluxspep = retr_numbspep(gdat.spectype)
    gdat.indxspep = [arange(gdat.numbspep[l]) for l in gdat.indxpopl]
    ## plotting
    ### number of bins for histogram plots of spectral parameters
    gdat.numbspepbins = 20
    ### number of standard deviations away from the mean of spectral parameters to plot
    gdat.numbstdvspepdist = 3.
    ### minima and maxima for spectral parameters
    gdat.numbstdv = 3.

    # number of components
    gdat.numbcomp = 3 + zeros(gdat.numbpopl, dtype=int)
    gdat.numbcompcolr = 3 + zeros(gdat.numbpopl, dtype=int)
    if gdat.numbener > 1:
        gdat.numbcomp += gdat.numbspep + gdat.numbener - 1
        gdat.numbcompcolr += gdat.numbspep
    gdat.maxmnumbcompcolr = amax(gdat.numbcompcolr)

    gdat.indxcomp = []
    for l in gdat.indxpopl:
        gdat.indxcomp.append(arange(gdat.numbcomp[l]))

    # total maximum number of PS
    gdat.maxmnumbpntstotl = sum(gdat.maxmnumbpnts)
    gdat.indxpntstotl = arange(gdat.maxmnumbpntstotl)
    gdat.maxmnumbpntscumr = cumsum(gdat.maxmnumbpnts)
    gdat.maxmnumbpntscuml = concatenate((array([0]), gdat.maxmnumbpntscumr[:-1]))
   
    # minimum number of PS
    gdat.minmnumbpnts = zeros(gdat.numbpopl, dtype=int)
   
    # maximum number of components
    gdat.maxmnumbcomp = gdat.maxmnumbpnts * gdat.numbcomp
    gdat.maxmnumbcompcumr = cumsum(gdat.maxmnumbcomp)
    gdat.maxmnumbcompcuml = concatenate((array([0]), gdat.maxmnumbcompcumr[:-1]))
    gdat.maxmnumbcomptotl = sum(gdat.maxmnumbcomp)
    
    # sweeps to be saved
    gdat.boolsave = zeros(gdat.numbswep, dtype=bool)
    gdat.indxswepsave = arange(gdat.numbburn, gdat.numbburn + gdat.numbsamp * gdat.factthin, gdat.factthin)
    gdat.boolsave[gdat.indxswepsave] = True
    gdat.indxsampsave = zeros(gdat.numbswep, dtype=int) - 1
    gdat.indxsampsave[gdat.indxswepsave] = arange(gdat.numbsamp)
    
    # log-prior register
    ## indices of penalization term
    indxlpripena = 0
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    if gdat.maxmnumbpntstotl > 0:
        gdat.numblpri = 2 + 2 * gdat.numbpopl
        if gdat.numbener > 1:
            gdat.numblpri += 2 * gdat.numbpopl
    else:
        gdat.numblpri = 0

    # set model sample vector indices
    retr_indxsamp(gdat, gdat.psfntype, gdat.spectype, gdat.varioaxi)

    if gdat.datatype == 'mock':
        # set mock sample vector indices
        retr_indxsamp(gdat, gdat.truepsfntype, gdat.truespectype, gdat.truevarioaxi, strgpara='true')

    # process index
    gdat.indxproc = arange(gdat.numbproc)

    # flag to indicate whether information from a deterministic catalog will be used or not
    # temp -- if datatype == 'inpt' trueinfo should depend on whether truexxxx are provided
    gdat.trueinfo = gdat.datatype == 'mock' or gdat.exprinfo
    
    # half size of the image where the sample catalog is compared against the reference
    gdat.maxmgangcomp = gdat.maxmgang * gdat.margfactcomp
    # half size of the spatial prior
    gdat.maxmgangmodl = gdat.maxmgang * gdat.margfactmodl

    # axes
    # temp
    gdat.liststrgpntspara = ['lgal', 'bgal'] + list(set([strg for strg in gdat.liststrgfluxspep[l] for l in gdat.indxpopl]))
    for strgpntspara in gdat.liststrgpntspara:
        setattr(gdat, 'numb' + strgpntspara + 'plot', 20)
   
    if gdat.pntstype == 'lens':
        gdat.minmmass = retr_massfrombein(gdat.minmflux)
        gdat.maxmmass = retr_massfrombein(gdat.maxmflux)
        retr_axis(gdat, 'bein', gdat.minmflux, gdat.maxmflux, 10)

    gdat.indxspepsind = 0
    gdat.indxspepcurv = 1
    gdat.indxspepexpo = 2

    gdat.numbsinddistpara = 2
    gdat.numbfluxdistpara = 4
    
    gdat.listlablcompfrac = ['Data']
    if gdat.pntstype == 'lght' or gdat.numbback > 1:
        gdat.listlablcompfrac.append('Total Model')
    if gdat.pntstype == 'lght':
        gdat.listlablcompfrac.append('PS')
    gdat.listlablcompfrac += gdat.lablback
    gdat.numblablcompfrac = len(gdat.listlablcompfrac)

    if gdat.strgfluxunit == None:
        gdat.strgfluxunitextn = ''
    else:
        gdat.strgfluxunitextn = ' [%s]' % gdat.strgfluxunit

    if gdat.numbener > 1:
        gdat.enerfluxdist = gdat.meanener[gdat.indxenerfluxdist]
        if gdat.enerfluxdist == 0.:
            raise Exception('Pivot energy cannot be zero.')
        gdat.enernorm = gdat.meanener / gdat.enerfluxdist
        gdat.factlogtenerpivt = log(gdat.enernorm)
        gdat.factspecener = gdat.enernorm**(-sqrt(amin(gdat.minmsinddistmean) * amax(gdat.maxmsinddistmean)))
        gdat.enerexpofact = gdat.enerfluxdist - gdat.meanener
    else:
        gdat.factspecener = array([1.])

    # angular deviation
    # temp -- check that gdat.numbangl does not degrade the performance
    if gdat.pntstype == 'lght':
        gdat.numbangl = 100
        retr_axis(gdat, 'angl', 0., gdat.maxmangl, gdat.numbangl)
        gdat.binsanglcosi = sort(cos(gdat.binsangl))
   
    gdat.meshbackener = meshgrid(gdat.indxback, gdat.indxener, indexing='ij')

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
    if gdat.strganglunit != '':
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
    
    # off-axis angle
    if gdat.varioaxi or gdat.truevarioaxi:
        gdat.numboaxi = 10
        gdat.minmoaxi = 0.
        gdat.maxmoaxi = 1.1 * sqrt(2.) * gdat.maxmgangmodl
        retr_axis(gdat, 'oaxi', gdat.minmoaxi, gdat.maxmoaxi, gdat.numboaxi)
        gdat.binsoaxiopen = gdat.binsoaxi[:-1]
    else:
        gdat.binsoaxi = None
        gdat.numboaxi = 1
    gdat.indxoaxi = arange(gdat.numboaxi)

    gdat.numbenerevtt = gdat.numbener * gdat.numbevtt
    
    # pixelization
    if gdat.datatype == 'mock':
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2

    gdat.numbchrototl = 5
    gdat.numbchrollik = 12

    # pivot off-axis scale
    gdat.oaxipivt = gdat.maxmgang

    # temp
    gdat.boolintpanglcosi = False

    # number of bins
    gdat.numbbins = 10

    # the function to measure time
    if gdat.strgfunctime == 'clck':
        gdat.functime = time.clock
    if gdat.strgfunctime == 'time':
        gdat.functime = time.time

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

    if gdat.pntstype == 'lens': 
        retr_axis(gdat, 'defl', -gdat.maxmgang, gdat.maxmgang, 50)

    # convenience variables
    gdat.indxfluxprox = arange(gdat.numbfluxprox)
    gdat.binsfluxprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
    gdat.meanfluxprox = sqrt(gdat.binsfluxprox[1:] * gdat.binsfluxprox[:-1])
    
    # temp
    if False:
        minm = None
        maxm = None
        for strgpara in ['sind']:
            for strgdata in ['', 'mock']:
                strgpopl = strgdata + 'numbpopl'
                numbpopl = getattr(gdat, strgpopl)
                for l in range(numbpopl):
                    if strgdata == '':
                        strgmean = strgdata + 'mean' + strgpara + 'distmean'
                        strgstdv = strgdata + 'stdv' + strgpara + 'distmean'
                    else:
                        strgmean = strgdata + strgpara + 'distmean'
                        strgstdv = strgdata + strgpara + 'distmean'
                    minm = min(minm, getattr(gdat, strgmean)[l] - gdat.numbstdv * getattr(gdat, strgstdv)[l])
                    maxm = max(maxm, getattr(gdat, strgmean)[l] + gdat.numbstdv * getattr(gdat, strgstdv)[l])
            setattr(gdat, 'minm' + strgpara, minm)
            setattr(gdat, 'maxm' + strgpara, maxm)
    
    else:
        gdat.minmsind = 1.
        gdat.maxmsind = 3.

    # temp
    gdat.minmcurv = 0.
    gdat.maxmcurv = 1.

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
    gdat.pathinpt = gdat.pathdata + 'inpt/'
    if gdat.datatype == 'inpt':
        
        path = gdat.pathinpt + gdat.strgexprflux
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
            gdat.sizepixl = sqrt(gdat.apix)
            gdat.indxpixlrofi = arange(gdat.numbsidecart**2)
            gdat.indxsidecart = arange(gdat.numbsidecart)
            gdat.indxsidemesh = meshgrid(gdat.indxsidecart, gdat.indxsidecart, indexing='ij')
            gdat.bgalgrid = gdat.bgalcart[gdat.indxsidemesh[1].flatten()]
            gdat.lgalgrid = gdat.lgalcart[gdat.indxsidemesh[0].flatten()]
            gdat.shapcart = (gdat.numbsidecart, gdat.numbsidecart)
            gdat.lgalgridcart = gdat.lgalgrid.reshape(gdat.shapcart)
            gdat.bgalgridcart = gdat.bgalgrid.reshape(gdat.shapcart)
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
    
    gdat.jcbnsplt = 2.**(2 - gdat.numbener)
    
    # power spectra
    if gdat.pixltype == 'cart':
        gdat.factkpcs = 1e-6
        gdat.numbwvecodim = 40
        gdat.indxwvecodim = arange(gdat.numbwvecodim)
        gdat.minmwvecodim = 2. * pi / sqrt(2) / gdat.maxmgang * gdat.factkpcs # [1/kpc]
        gdat.maxmwvecodim = 2. * pi / gdat.sizepixl * gdat.factkpcs
        gdat.binswvecodim = linspace(gdat.minmwvecodim, gdat.maxmwvecodim, gdat.numbwvecodim + 1)
        gdat.meanwvecodim = (gdat.binswvecodim[1:] + gdat.binswvecodim[:-1]) / 2.
        gdat.numbsidewvec = gdat.numbsidecart / 2
        temp = fft.fftfreq(gdat.numbsidewvec, gdat.sizepixl)
        gdat.meanwveclgal, gdat.meanwvecbgal = meshgrid(temp, temp, indexing='ij')
        gdat.meanwveclgal *= gdat.factkpcs
        gdat.meanwvecbgal *= gdat.factkpcs
        gdat.meanwvec = sqrt(gdat.meanwveclgal**2 + gdat.meanwvecbgal**2)

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
            path = gdat.pathinpt + gdat.strgexpo
            gdat.expo = pf.getdata(path)
            if amin(gdat.expo) == amax(gdat.expo):
                raise Exception('Bad input exposure map.')
                return
            if gdat.pixltype == 'cart':
                gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))

    # backgrounds
    gdat.backflux = []
    for c in gdat.indxback:
        if isinstance(gdat.back[c], float):
            if gdat.datatype == 'mock':
                if gdat.pixltype == 'heal':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + gdat.back[c]
                if gdat.pixltype == 'cart':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + gdat.back[c]
                if gdat.pixltype == 'unbd':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbdatasamp, gdat.numbevttfull)) + gdat.back[c]
            if gdat.datatype == 'inpt':
                backfluxtemp = zeros_like(gdat.exprdataflux) + gdat.back[c]
        else:
            path = gdat.pathinpt + gdat.back[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'cart':
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        gdat.backflux.append(backfluxtemp)
  
    if gdat.pixltype == 'cart':
        gdat.backfluxcart = []
        for c in gdat.indxback:
            gdat.backfluxcart.append(gdat.backflux[c].reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt)))
        gdat.expocart = gdat.expo.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))

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
        path = gdat.pathpixlcnvt + 'pixlcnvt_%09g.p' % gdat.maxmgang

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

    # temp
    if gdat.pntstype == 'lens':
        print 'gdat.numbsidecart'
        print gdat.numbsidecart
        print 'gdat.backflux[0]'
        print gdat.backflux[0].shape

        gdat.backfluxlens = gdat.backflux[0][0, :, 0].reshape((gdat.numbsidecart, gdat.numbsidecart))
    
    if gdat.correxpo:
        gdat.backcnts = []
        gdat.backcntstotl = zeros_like(gdat.expo)
        for c in gdat.indxback:
            backcntstemp = retr_cntsmaps(gdat, gdat.backflux[c])
            gdat.backcnts.append(backcntstemp)
            gdat.backcntstotl[:] += backcntstemp 
    
    if False and gdat.evalcirc and gdat.numbpixl * gdat.maxmnumbpntstotl < 1e5:
        gdat.calcerrr = True
    else:
        gdat.calcerrr = False
    
    if gdat.pntstype == 'lght':
        gdat.exprpsfn = retr_psfn(gdat, gdat.exprpsfp, gdat.indxener, gdat.binsangl, gdat.exprpsfntype, gdat.binsoaxi, gdat.exprvarioaxi)

    if gdat.evalcirc:
        # determine the maximum angle at which the PS flux map will be computed
        gdat.maxmangleval = empty(gdat.numbfluxprox)
        for h in gdat.indxfluxprox:
            if gdat.specfraceval == 0:
                gdat.maxmangleval[h] = 3. * gdat.maxmgangmodl
            else:   
                frac = gdat.specfraceval * gdat.binsfluxprox[0] / gdat.binsfluxprox[h+1]
                psfnwdth = retr_psfnwdth(gdat, gdat.exprpsfn, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathpixlprox + 'indxpixlprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), \
                                                                                                            1e2 * amax(gdat.maxmangleval), gdat.numbfluxprox)
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
                dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
                dist[j] = 0.
                for h in range(gdat.numbfluxprox):
                    indxpixlproxtemp = where(dist < gdat.maxmangleval[h])[0]
                    gdat.indxpixlprox[h].append(indxpixlproxtemp)
                cntrsave = tdpy.util.show_prog(j, gdat.indxpixl.size, cntrsave)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()

    # plot settings
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
    ## marker opacity
    gdat.alphmrkr = 0.5
    gdat.alphpnts = 0.4
    gdat.alphmaps = 1.
    
    gdat.numbbinsplot = 10

    ## flux
    gdat.minmfluxplot = gdat.minmflux
    gdat.maxmfluxplot = gdat.maxmflux
    if gdat.trueinfo:
        if gdat.trueminmflux != None:
            gdat.minmfluxplot = min(gdat.minmfluxplot, gdat.trueminmflux)
        if gdat.truemaxmflux != None:
            gdat.maxmfluxplot = max(gdat.maxmfluxplot, gdat.truemaxmflux)
    retr_axis(gdat, 'fluxplot', gdat.minmfluxplot, gdat.maxmfluxplot, gdat.numbbinsplot, scal='logt')
    
    # plotting
    gdat.numbtickcbar = 11
    gdat.minmconv = 1e-2
    gdat.maxmconv = 1e1
    gdat.minmdeflcomp = 1.
    gdat.maxmdeflcomp = 1.1

    if gdat.pntstype == 'lens':
        liststrgcbar = ['conv', 'deflcomp']
        for strgcbar in liststrgcbar:
            retr_ticklabl(gdat, strgcbar)
    

def setpfinl(gdat, boolinitsetp=False):

    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = retr_cntsmaps(gdat, gdat.exprdataflux)
    
    if gdat.pixltype != 'unbd':
        gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix
        if gdat.enerbins:
            gdat.datafluxmean /= gdat.diffener
    else:
        gdat.datafluxmean = array([gdat.numbdatasamp / gdat.apix])
        
    if gdat.pntstype == 'lght':
        gdat.limsangl = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        gdat.limspsfn = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.truevarioaxi:
                    psfn = gdat.exprpsfn[i, :, m, 0]
                else:
                    psfn = gdat.exprpsfn[i, :, m]
                maxmpsfn = amax(psfn)
                gdat.limsangl[i][m] = [0., gdat.binsangl[amax(where(psfn > 1e-6 * maxmpsfn)[0])] * gdat.anglfact]
                gdat.limspsfn[i][m] = [maxmpsfn * 1e-6, maxmpsfn]
       
    # get the experimental catalog
    if gdat.exprinfo:
        
        gdat.exprcnts = None
        if gdat.exprtype == 'ferm':
            retr_fermdata(gdat)
        if gdat.exprtype == 'chan':
            retr_chandata(gdat)
    
        # rotate PS coordinates to the ROI center
        if gdat.lgalcntr != 0. or gdat.bgalcntr != 0.:
            rttr = hp.rotator.Rotator(rot=[rad2deg(gdat.lgalcntr), rad2deg(gdat.bgalcntr), 0.], deg=True, eulertype='ZYX')
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
    
        # reorder PS with respect to flux
        indxpnts = argsort(gdat.exprspec[0, gdat.indxenerfluxdist[0], :])[::-1]
        gdat.exprlgal = gdat.exprlgal[indxpnts]
        gdat.exprbgal = gdat.exprbgal[indxpnts]
        gdat.exprspec[0, :, :] = gdat.exprspec[0, :, indxpnts].T
        if gdat.exprcnts != None:
            gdat.exprcnts = gdat.exprcnts[:, indxpnts, :]

        # compute the catalog counts based on the exposure
        gdat.exprcntscalc = empty((gdat.numbener, gdat.exprnumbpnts, gdat.numbevtt))
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                indxpixltemp = retr_indxpixl(gdat, gdat.exprbgal, gdat.exprlgal)
                gdat.exprcntscalc[i, :, m] = gdat.exprspec[0, i, :] * gdat.expo[i, indxpixltemp, m] * gdat.diffener[i]
      
        if gdat.strgcnfg == 'pcat_chan_inpt':
            print 'gdat.exprcnts'
            print gdat.exprcnts[:, :, 0].T
            print 'gdat.exprcntscalc'
            print gdat.exprcntscalc[:, :, 0].T

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

    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(10000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)
    gdat.indxcubesave = meshgrid(gdat.indxener, gdat.indxpixlsave, gdat.indxevtt, indexing='ij')
    
    if gdat.correxpo:
        # limits on counts, which are used to bin or overplot PS counts 
        # temp
        gdat.minmcnts = 0.1 * gdat.minmflux * mean(mean(gdat.expo, 1), 1)
        gdat.maxmcnts = gdat.maxmflux * mean(mean(gdat.expo, 1), 1)
        if gdat.enerbins:
            gdat.minmcnts *= gdat.diffener * gdat.factspecener
            gdat.maxmcnts *= gdat.diffener * gdat.factspecener
        gdat.binscnts = zeros((gdat.numbener, gdat.numbbinsplot + 1))
        for i in gdat.indxener:
            gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbbinsplot + 1) # [1]
        gdat.meancnts = sqrt(gdat.binscnts[:, :-1] * gdat.binscnts[:, 1:]) 

    # load the true data into the reference data structure
    if gdat.trueinfo:
        ## mock data
        if gdat.datatype == 'mock':
            for k in gdat.indxfixp:
                
                # temp -- allow mismodeling
                continue

                if gdat.namefixp[k][:-1].endswith('pop'):
                    strgpopl = gdat.namefixp[k][-1]
                    l = int(strgpopl)
                    #if gdat.truenumbpopl != gdat.numbpopl or gdat.namefixp[k].startswith('fluxdist') and gdat.truefluxdisttype[l] != gdat.fluxdisttype:
                    #    continue
                    strg = gdat.namefixp[k][:-4]
                    print 'true' + strg
                    gdat.corrfixp[k] = getattr(gdat, 'mock' + strg)[l]
                else:
                    strg = gdat.namefixp[k]
                    print 'true' + strg
                    gdat.corrfixp[k] = getattr(gdat, 'mock' + strg)
    
        ## Real data
        if gdat.datatype == 'inpt':
            gdat.truefixp[gdat.indxfixpnumbpnts] = array([gdat.exprnumbpnts], dtype=int)
            gdat.truefixp[gdat.indxfixpmeanpnts] = gdat.truefixp[gdat.indxfixpnumbpnts]
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
            
    if gdat.pixltype == 'unbd':
        gdat.bgalgrid = gdat.datacnts[0, :, 0, 0]
        gdat.lgalgrid = gdat.datacnts[0, :, 0, 1]
   
    # temp -- needed for experimental PSF
    # if gdat.evalpsfnpnts:
    #     gdat.truenumbpsfpform, gdat.truenumbpsfpoaxi, gdat.truenumbpsfptotl, gdat.trueindxpsfpoaxinorm, gdat.trueindxpsfpoaxiindx = \
    #                                                                                         retr_indxpsfp(gdat, gdat.truepsfntype, gdat.truevarioaxi)
    #     if gdat.truevarioaxi:
    #         gdat.truefactoaxi = retr_factoaxi(gdat, gdat.binsoaxi, gdat.truepsfp[gdat.trueindxpsfpoaxinorm], gdat.truepsfp[gdat.trueindxpsfpoaxiindx])
    
    # get count data
    if gdat.pixltype == 'cart':
        # temp
        gdat.indxxaximaxm, gdat.indxyaximaxm = tdpy.util.retr_indximagmaxm(gdat.datacnts[0, :, 0].reshape((gdat.numbsidecart, gdat.numbsidecart)))

    if gdat.verbtype > 1 and boolinitsetp:
        print 'fixp'
        if gdat.datatype == 'mock':
            liststrgpara = ['', 'true']
        else:
            liststrgpara = ['']
        for strgpara in liststrgpara:
            if strgpara == '':
                print 'modl'
            else:
                print strgpara
            if strgpara == '':
                listfeat = ('name', 'strg', 'scal', 'minm', 'maxm', 'true')
                print '%20s%20s%5s%20s%20s%20s' % listfeat
            else:
                listfeat = ('name', 'strg', 'scal', 'minm', 'maxm', 'sampvarb')
                print '%20s%20s%5s%20s%20s%20s' % listfeat
            for k in getattr(gdat, strgpara + 'indxfixp'):
                if strgpara == '':
                    if gdat.truefixp[k] == None:
                        strg = '%20s' % 'None'
                    else:
                        strg = '%20.6g' % gdat.truefixp[k]
                    print '%20s%20s%5s%20.6g%20.6g%s' % (gdat.namefixp[k], gdat.strgfixp[k], gdat.scalfixp[k], gdat.minmfixp[k], gdat.maxmfixp[k], strg)
                else:
                    print '%20s%20s%5s%20.6g%20.6g%20.6g' % (gdat.namefixp[k], gdat.strgfixp[k], gdat.scalfixp[k], gdat.minmfixp[k], gdat.maxmfixp[k], gdat.truefixp[k])
        
    if gdat.trueinfo and gdat.correxpo and gdat.pntstype == 'lght':
        truebackcnts = []
        gdat.truesigm = []
        if gdat.numbtrap > 0:
            for l in gdat.indxpopl:
                indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
                truebackcntstemp = zeros((gdat.numbener, gdat.truefixp[gdat.indxfixpnumbpnts[l]], gdat.numbevtt))
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

        if gdat.numbtrap > 0:
            gdat.truefluxbrgt, gdat.truefluxbrgtassc = retr_fluxbrgt(gdat, concatenate(gdat.truelgal), concatenate(gdat.truebgal), \
                                                                                                        concatenate(gdat.truespec, axis=2)[0, gdat.indxenerfluxdist[0], :])
    
   
    # proposals

    
    liststrg = []
    if gdat.numbtrap > 0:
        liststrg += ['numbpnts', 'meanpnts']
        for l in gdat.indxpopl:
            if gdat.fluxdisttype == 'powr':  
                liststrg += ['fluxdistslop']
            if gdat.fluxdisttype == 'brok':  
                liststrg += ['fluxdistbrek', 'fluxdistsloplowr', 'fluxdistslopuppr']
            if gdat.numbener > 1:
                liststrg += ['sinddistmean', 'sinddiststdv']
    liststrg += ['psfp', 'bacp']
    if gdat.pntstype == 'lens':
        liststrg += ['lenp']
        
    gdat.indxfixpactv = array([], dtype=int)
    gdat.indxfixpactvprop = array([], dtype=int)
    for strgtype in ['indxfixpactv', 'indxfixpactvprop']:
        for strg in liststrg:
            
            if strgtype == 'indxfixpactvprop':
                if strg[4:8] == 'dist' or strg == 'meanpnts':
                    strgtemp = 'hypr'
                else:
                    strgtemp = strg
                thisbool = getattr(gdat, 'prop' + strgtemp) and strg != 'numbpnts'
            else:
                thisbool = True
            if thisbool:
                setattr(gdat, strgtype, append(getattr(gdat, strgtype), getattr(gdat, 'indxfixp' + strg)))
        
    
    gdat.numbfixpactvprop = gdat.indxfixpactvprop.size
    gdat.indxactvprop = arange(gdat.numbfixpactvprop)
    
    gdat.indxfixpdist = setdiff1d(gdat.indxfixphypr, gdat.indxfixpmeanpnts)
    gdat.indxfixpdistactv = intersect1d(gdat.indxfixpdist, gdat.indxfixpactv) 
    gdat.indxfixphypractv = intersect1d(gdat.indxfixphypr, gdat.indxfixpactv)

    gdat.indxfixpiact = setdiff1d(gdat.indxfixp, gdat.indxfixpactvprop)
    gdat.numbfixpiact = gdat.indxfixpiact.size
    gdat.indxiact = arange(gdat.numbfixpiact)
    gdat.indxsampiact = zeros(gdat.numbfixpiact)
    cntr = 1
    for k in gdat.indxfixp:
        if k in gdat.indxsampiact:
            gdat.indxsampiact[cntr] = k
            cntr += 1

    gdat.strgprop = array([])
    gdat.nameprop = array([])
   
    if gdat.probtran == None:
        if gdat.numbtrap > 0:
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
       
    cntr = tdpy.util.cntr()
    if gdat.numbtrap > 0.:
    
        if gdat.probtran > 0.:
            # birth
            gdat.indxpropbrth = cntr.incr()
            gdat.strgprop = append(gdat.strgprop, r'$\mathcal{B}$')
            gdat.nameprop = append(gdat.nameprop, 'brth')
            
            # death
            gdat.indxpropdeth = cntr.incr()
            gdat.strgprop = append(gdat.strgprop, r'$\mathcal{D}$')
            gdat.nameprop = append(gdat.nameprop, 'deth')
            
            if gdat.probbrde < 1.:
                # split
                gdat.strgprop = append(gdat.strgprop, r'$\mathcal{S}$')
                gdat.nameprop = append(gdat.nameprop, 'splt')
                gdat.indxpropsplt = cntr.incr()
                
                # merge
                gdat.strgprop = append(gdat.strgprop, r'$\mathcal{M}$')
                gdat.nameprop = append(gdat.nameprop, 'merg')
                gdat.indxpropmerg = cntr.incr()
    
    gdat.numbprop = gdat.strgprop.size
    gdat.indxprop = arange(gdat.numbprop)
   
    gdat.indxactvconv = zeros(gdat.numbfixp, dtype=int)
    gdat.indxactvconv[gdat.indxfixpactvprop] = arange(gdat.numbfixpactvprop, dtype=int)
    
    gdat.strgcompcolr = ['lgal', 'bgal', 'flux']
    gdat.indxstdplgal = gdat.numbfixpactvprop
    gdat.indxstdpbgal = gdat.numbfixpactvprop + 1
    gdat.indxstdpflux = gdat.numbfixpactvprop + 2
    if gdat.numbener > 1:
        gdat.strgcompcolr += ['spep']
        gdat.indxstdpspep = gdat.numbfixpactvprop + 3
    gdat.numbstdp = gdat.numbfixpactvprop + gdat.maxmnumbcompcolr
    gdat.strgstdp = concatenate((gdat.strgfixp[gdat.indxfixpactvprop], gdat.strgcompcolr))
    gdat.namestdp = concatenate((gdat.namefixp[gdat.indxfixpactvprop], gdat.strgcompcolr))
    gdat.indxstdp = arange(gdat.numbstdp)
    gdat.indxstdpfixp = arange(gdat.numbfixpactvprop)
    gdat.indxstdpcomp = setdiff1d(gdat.indxstdp, gdat.indxstdpfixp)

   
    gdat.indxstdpactv = gdat.indxfixpactvprop
    if gdat.numbtrap > 0:
        gdat.indxstdpactv = concatenate((gdat.indxstdpactv, array([gdat.indxstdplgal, gdat.indxstdpbgal, gdat.indxstdpflux])))
        if gdat.numbener > 1:
            gdat.indxstdpactv = concatenate((gdat.indxstdpactv, array([gdat.indxstdpspep])))

    # proposal scale
    gdat.stdvstdp = 1e-4 + zeros(gdat.numbstdp)
    
    listindxsampunsd = []
    numbsampcumu = 0
    for l in gdat.indxpopl:
        for k in gdat.indxpnts[l]:
            listindxsampunsd.append(gdat.indxsampcomp[0] + numbsampcumu + k * gdat.numbcomp[l] + gdat.indxcompunsd[l])
        numbsampcumu += gdat.maxmnumbpnts[l] * gdat.numbcomp[l]
    if gdat.numbtrap > 0:
        gdat.indxsampunsd = concatenate(listindxsampunsd)
    else:
        gdat.indxsampunsd = []

    ## color
    gdat.binssind = linspace(gdat.minmsind, gdat.maxmsind, gdat.numbspepbins + 1)
    gdat.meansind = (gdat.binssind[1:] + gdat.binssind[:-1]) / 2.
    gdat.diffsind = gdat.binssind[1:] - gdat.binssind[:-1]

    gdat.minmspecplot = gdat.minmfluxplot * gdat.factspecener
    gdat.maxmspecplot = gdat.maxmfluxplot * gdat.factspecener
    gdat.binsspecplot = gdat.binsfluxplot[None, :] * gdat.factspecener[:, None]
    gdat.meanspecplot = empty((gdat.numbener, gdat.numbbinsplot))
    for i in gdat.indxener:
        gdat.meanspecplot[i, :] = sqrt(gdat.binsspecplot[i, 1:] * gdat.binsspecplot[i, :-1])

    # determine the indices of true point sources, which will be compared againts the model sources
    if gdat.trueinfo:
        gdat.trueindxpntscomp = []
        for l in gdat.indxpopl:
            trueindxpntstemp = where((fabs(gdat.truelgal[l]) < gdat.maxmgangcomp) & (fabs(gdat.truebgal[l]) < gdat.maxmgangcomp))[0]
            gdat.trueindxpntscomp.append(trueindxpntstemp)

    # sanity checks
    # temp
    if (fabs(gdat.datacnts - rint(gdat.datacnts)) > 1e-3).any() and boolinitsetp:
        print 'Fractional counts!'

    if amin(gdat.datacnts) < 0. and boolinitsetp:
        print 'Negative counts!'

    retr_datatick(gdat)

    if gdat.verbtype > 1 and boolinitsetp:
        if gdat.pntstype == 'lght' and gdat.pixltype != 'unbd':
            print 'Memory budget: indxpixlprox'
            totl = 0.
            for h in gdat.indxfluxprox:
                for n in gdat.indxpixl:
                    totl += sys.getsizeof(gdat.indxpixlprox[h][n]) / 2.**20
            print '%.4g MB' % totl
    

def retr_datatick(gdat):

    # data count limits
    gdat.minmdatacnts = 0.
    gdat.maxmdatacnts = amax(sum(gdat.datacnts, 2))
    gdat.maxmresicnts = ceil(gdat.maxmdatacnts * 0.1)
    gdat.minmresicnts = -gdat.maxmresicnts
    
    # plotting
    liststrgcbar = ['datacnts', 'resicnts']
    for strgcbar in liststrgcbar:
        retr_ticklabl(gdat, strgcbar)


def retr_ticklabl(gdat, strgcbar):

    if gdat.scalmaps == 'asnh':
        for strglimt in ['minm', 'maxm']:
            setattr(gdat, strglimt + strgcbar, arcsinh(getattr(gdat, strglimt + strgcbar)))

    tick = linspace(getattr(gdat, 'minm' + strgcbar), getattr(gdat, 'maxm' + strgcbar), gdat.numbtickcbar)
    
    labl = empty(gdat.numbtickcbar, dtype=object)
    for k in range(gdat.numbtickcbar):
        if gdat.scalmaps == 'asnh':
            labl[k] = '%.3g' % sinh(tick[k])
        else:
            labl[k] = '%.3g' % tick[k]
    
    setattr(gdat, 'tick' + strgcbar, tick)
    setattr(gdat, 'labl' + strgcbar, labl)
    

def retr_fromgdat(gdat, gdatmodi, strg, strgvarb, errr=False):
    
    if strg == 'this':
        varb = getattr(gdatmodi, strg + strgvarb)
    elif strg == 'post':
        if errr:
            varb = getattr(gdat, 'errr' + strgvarb)
        else:
            varb = getattr(gdat, 'medi' + strgvarb)
    else:
        varb = getattr(gdat, strg + strgvarb)

    return varb


def retr_indxsamp(gdat, psfntype, spectype, varioaxi, strgpara=''):

    numbback = getattr(gdat, strgpara + 'numbback')
    indxback = arange(numbback)
    
    numbbacp = numbback * gdat.numbener
    
    # population index vector
    numbpopl = getattr(gdat, strgpara + 'numbpopl')
    indxpopl = arange(numbpopl, dtype=int) 
    
    cntr = tdpy.util.cntr()
    
    indxfixpfluxdistbrek = []
    if gdat.maxmnumbpntstotl > 0:
        indxfixpnumbpnts = arange(numbpopl) + cntr.incr(numbpopl)
        indxfixpmeanpnts = arange(numbpopl) + cntr.incr(numbpopl)
        indxfixphypr = [indxfixpnumbpnts, indxfixpmeanpnts] 
        if gdat.fluxdisttype == 'powr':
            indxfixpfluxdistslop = arange(numbpopl) + cntr.incr(numbpopl)
            indxfixphypr += [indxfixpfluxdistslop]
        else:
            indxfixpfluxdistbrek = arange(numbpopl) + cntr.incr(numbpopl)
            indxfixpfluxdistsloplowr = arange(numbpopl) + cntr.incr(numbpopl)
            indxfixpfluxdistslopuppr = arange(numbpopl) + cntr.incr(numbpopl)
            indxfixphypr += [indxfixpfluxdistbrek, indxfixpfluxdistsloplowr, indxfixpfluxdistslopuppr]
        if gdat.numbener > 1:
            indxfixpsinddistmean = arange(numbpopl) + cntr.incr(numbpopl)
            indxfixpsinddiststdv = arange(numbpopl) + cntr.incr(numbpopl)
            indxfixphypr += [indxfixpsinddistmean, indxfixpsinddiststdv]
            if spectype == 'curv':
                indxfixpcurvdistmean = arange(numbpopl) + cntr.incr(numbpopl)
                indxfixpcurvdiststdv = arange(numbpopl) + cntr.incr(numbpopl)
                indxfixphypr += [indxfixpcurvdistmean, indxfixpcurvdiststdv]
            if spectype == 'expo':
                indxfixpexpodistmean = arange(numbpopl) + cntr.incr(numbpopl)
                indxfixpexpodiststdv = arange(numbpopl) + cntr.incr(numbpopl)
                indxfixphypr += [indxfixpexpodistmean, indxfixpexpodiststdv]
             
        indxfixphypr = concatenate(indxfixphypr)

    indxfixpsigc = []
    indxfixpsigt = []
    indxfixpgamc = []
    indxfixpgamt = []
    indxfixppsff = []
    indxfixpoaxinorm = []
    indxfixpoaxiindx = []
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if psfntype == 'singgaus':
                indxfixpsigc.append(cntr.incr())
            if psfntype == 'singking':
                indxfixpsigc.append(cntr.incr())
                indxfixpgamc.append(cntr.incr())
            if psfntype == 'doubgaus':
                indxfixpsigc.append(cntr.incr())
                indxfixpsigt.append(cntr.incr())
                indxfixppsff.append(cntr.incr())
            if psfntype == 'gausking':
                indxfixpsigc.append(cntr.incr())
                indxfixpgamc.append(cntr.incr())
                indxfixpsigt.append(cntr.incr())
                indxfixppsff.append(cntr.incr())
            if psfntype == 'doubking':
                indxfixpsigc.append(cntr.incr())
                indxfixpgamc.append(cntr.incr())
                indxfixpsigt.append(cntr.incr())
                indxfixpgamt.append(cntr.incr())
                indxfixppsff.append(cntr.incr())
            if psfntype == 'neww':
                pass
            if varioaxi:
                indxfixpoaxinorm.append(cntr.incr())
                indxfixpoaxiindx.append(cntr.incr())
    
    indxfixpsigm = indxfixpsigc + indxfixpsigt
    indxfixpgamm = indxfixpgamc + indxfixpgamt
    indxfixppsfp = indxfixpsigc + indxfixpsigt + indxfixpgamc + indxfixpgamt + indxfixppsff + indxfixpoaxinorm + indxfixpoaxiindx
    indxfixppsfp = sort(array(indxfixppsfp))

    numbpsfpform, numbpsfpoaxi, numbpsfptotl, indxpsfpoaxinorm, indxpsfpoaxiindx = retr_indxpsfp(gdat, psfntype, varioaxi)
    
    numbpsfptotlevtt = gdat.numbevtt * numbpsfptotl
    numbpsfptotlener = gdat.numbener * numbpsfptotl
    numbpsfp = numbpsfptotl * gdat.numbener * gdat.numbevtt
    indxpsfpoaxi = arange(numbpsfpoaxi) 
    indxpsfpform = arange(numbpsfpform)
    indxpsfptotl = arange(numbpsfptotl)
   
    if varioaxi:
        indxfixppsfpoaxinorm = indxfixppsfp[0] + indxpsfpoaxinorm
        indxfixppsfpoaxiindx = indxfixppsfp[0] + indxpsfpoaxiindx
        indxfixppsfpoaxi = sort(concatenate((indxfixppsfpoaxinorm, indxfixppsfpoaxiindx)))

    indxpsfp = arange(numbpsfp)
    indxpsfpinit = numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)

    indxfixpbacp = arange(numbbacp) + cntr.incr(numbbacp)

    indxfixplenp = []
    indxfixpanglsour = []
    indxfixpanglhost = []
    indxfixpangllens = []
    indxfixpsour = []
    indxfixphost = []
    if gdat.pntstype == 'lens':
        indxfixplgalsour = cntr.incr()
        indxfixpbgalsour = cntr.incr()
        indxfixpspecsour = cntr.incr()
        indxfixpsizesour = cntr.incr()
        indxfixpellpsour = cntr.incr()
        indxfixpanglsour = cntr.incr()
        indxfixplgalhost = cntr.incr()
        indxfixpbgalhost = cntr.incr()
        indxfixpspechost = cntr.incr()
        indxfixpsizehost = cntr.incr()
        indxfixpbeinhost = cntr.incr()
        indxfixpellphost = cntr.incr()
        indxfixpanglhost = cntr.incr()
        indxfixpsherhost = cntr.incr()
        indxfixpsanghost = cntr.incr()
        indxfixpsour = [indxfixplgalsour, indxfixpbgalsour, indxfixpspecsour, indxfixpsizesour, indxfixpellpsour, indxfixpanglsour]
        indxfixphost = [indxfixplgalhost, indxfixpbgalhost, indxfixpspechost, indxfixpsizehost, indxfixpbeinhost, indxfixpellphost, \
                                                                                                indxfixpanglhost, indxfixpsherhost, indxfixpsanghost]
        indxfixpemishost = [indxfixplgalhost, indxfixpbgalhost, indxfixpspechost, indxfixpellphost, indxfixpanglhost]
        indxfixplenp = list(set(indxfixpsour + indxfixphost + indxfixpemishost))
    
    # number of fixed-dimension parameters
    numbfixp = cntr.incr(0)
    # indices of fixed-dimension parameters
    indxfixp = arange(numbfixp)

    # total number of parameters
    numbpara = numbfixp + gdat.maxmnumbcomptotl
    indxsampcomp = arange(numbfixp, numbpara)
    indxpara = arange(numbpara)
    
    # transdimensional parameters
    numbtrap = numbpara - numbfixp

    # construct the fixed parameter structure
    namefixp = zeros(numbfixp, dtype=object)
    strgfixp = zeros(numbfixp, dtype=object)
    strgfixpunit = zeros(numbfixp, dtype=object)
    scalfixp = zeros(numbfixp, dtype=object)
    minmfixp = zeros(numbfixp)
    maxmfixp = zeros(numbfixp)
    factfixp = zeros(numbfixp)
    meanfixp = zeros(numbfixp)
    stdvfixp = zeros(numbfixp)
    cdfnminmfixp = empty(numbfixp)
    cdfndifffixp = empty(numbfixp)
    factfixpplot = ones(numbfixp)

    for k in indxfixp:
        if k in indxfixpnumbpnts or k in indxfixphypr:
            l = k % numbpopl

            if k in indxfixpnumbpnts:
                namefixp[k] = 'numbpntspop%d' % l
                strgfixp[k] = '$N$'
                scalfixp[k] = 'pois'
                
            if k in indxfixpmeanpnts:
                namefixp[k] = 'meanpntspop%d' % l
                strgfixp[k] = r'$\mu$'
                scalfixp[k] = 'logt'
    
            if gdat.fluxdisttype == 'powr':
                if k in indxfixpfluxdistslop:
                    namefixp[k] = 'fluxdistsloppop%d' % l
                    strgfixp[k] = r'$\alpha$'
                    scalfixp[k] = 'atan'
            else: 
                if k in indxfixpfluxdistbrek:
                    namefixp[k] = 'fluxdistbrekpop%d' % l
                    strgfixp[k] = '$f_b$'
                    scalfixp[k] = 'logt'
    
                if k in indxfixpfluxdistsloplowr:
                    namefixp[k] = 'fluxdistsloplowrpop%d' % l
                    strgfixp[k] = r'$\alpha_l$'
                    scalfixp[k] = 'atan'
    
                if k in indxfixpfluxdistslopuppr:
                    namefixp[k] = 'fluxdistslopupprpop%d' % l
                    strgfixp[k] = r'$\alpha_u$'
                    scalfixp[k] = 'atan'
            
            if gdat.numbener > 1:
                if k in indxfixpsinddistmean:
                    namefixp[k] = 'sinddistmeanpop%d' % l
                    strgfixp[k] = r'$\lambda_{s}$'
                    scalfixp[k] = 'atan'

                if k in indxfixpsinddiststdv:
                    namefixp[k] = 'sinddiststdvpop%d' % l
                    strgfixp[k] = r'$\sigma_{s}$'
                    scalfixp[k] = 'logt'

                if gdat.spectype == 'curv':
                    if k in indxfixpcurvdistmean:
                        namefixp[k] = 'curvdistmeanpop%d' % l
                        strgfixp[k] = r'$\lambda_{c}$'
                        scalfixp[k] = 'atan'

                    if k in indxfixpcurvdiststdv:
                        namefixp[k] = 'curvdiststdvpop%d' % l
                        strgfixp[k] = r'$\sigma_{c}$'
                        scalfixp[k] = 'logt'

                if gdat.spectype == 'curv':
                    if k in indxfixpexpodistmean:
                        namefixp[k] = 'expodistmeanpop%d' % l
                        strgfixp[k] = r'$\lambda_{\epsilon}$'
                        scalfixp[k] = 'logt'

                    if k in indxfixpexpodiststdv:
                        namefixp[k] = 'expodiststdvpop%d' % l
                        strgfixp[k] = r'$\sigma_{\epsilon}$'
                        scalfixp[k] = 'logt'

        if k in indxfixppsfp:
            if gdat.psfninfoprio:
                scalfixp[k] = 'gaus'
                n = k - indxfixppsfp[0]
                meanfixp[k] = gdat.meanpsfp[n]
                stdvfixp[k] = gdat.meanpsfp[n]
            else:
                if k in indxfixpsigm:
                    scalfixp[k] = 'logt'
                if k in indxfixpgamm:
                    scalfixp[k] = 'atan'
                if k in indxfixppsff:
                    scalfixp[k] = 'atan'
                if k in indxfixpoaxinorm:
                    scalfixp[k] = 'logt'
                if k in indxfixpoaxiindx:
                    scalfixp[k] = 'atan'
                
            # strings for PSF parameters
            if k in indxfixpsigm:
                strgvarbtemp = '\sigma'
                strgnametemp = 'sigm'
                factfixpplot[k] = gdat.anglfact
            if k in indxfixpgamm:
                strgvarbtemp = '\gamma'
                strgnametemp = 'gamm'
            if k in indxfixppsff:
                strgvarbtemp = 'f'
                strgnametemp = 'psff'
            if k in indxfixpoaxinorm:
                strgvarbtemp = 'a'
                strgnametemp = 'oaxinorm'
            if k in indxfixpoaxiindx:
                strgvarbtemp = 'b'
                strgnametemp = 'oaxiindx'
            if (k in indxfixpsigc or k in indxfixpgamc) and psfntype == 'doubgaus' or psfntype == 'gausking' or psfntype == 'doubking':
                    strgcomptemp = 'c'
                    strgnametemp = strgnametemp[:-1] + 'c'
            elif (k in indxfixpsigt or k in indxfixpgamt) and psfntype == 'gausking' or psfntype == 'doubking':
                    strgcomptemp = 't'
                    strgnametemp = strgnametemp[:-1] + 't'
            else:
                strgcomptemp = ''
            if gdat.numbener > 1:
                indxenertemp = gdat.indxenerincl[((k - indxfixppsfp[0]) % (gdat.numbener * numbpsfptotl)) // numbpsfptotl]
                strgenertemp = '%s' % indxenertemp
            else:
                strgenertemp = ''
            if gdat.numbevtt > 1:
                indxevtttemp = gdat.indxevttincl[(k - indxfixppsfp[0]) // (gdat.numbener * numbpsfptotl)]
                strgevtttemp = '%s' % indxevtttemp
            else:
                strgevtttemp = ''
            namefixp[k] = '%s%s%s' % (strgnametemp, strgenertemp, strgevtttemp)
            strgfixp[k] = r'$%s^{%s}_{%s%s}$' % (strgvarbtemp, strgcomptemp, strgenertemp, strgevtttemp)
        
        if k in indxfixpbacp:
            c = (k - indxfixpbacp[0]) // numbback
            if gdat.numbener > 1:
                i = (k - indxfixpbacp[0]) % numbback
                strgenertemp = '%d' % i
            else:
                strgenertemp = ''

            if numbback > 1:
                strgbacktemp = '%d' % c
            else:
                strgbacktemp = ''
            namefixp[k] = 'bacp'
            strgfixp[k] = '$A_{%s%s}$' % (strgbacktemp, strgenertemp)
            scalfixp[k] = 'logt'
        
        if gdat.pntstype == 'lens':
            if k in indxfixplenp:
                if k == indxfixplgalsour:
                    namefixp[k] = 'lgalsour'
                    strgfixp[k] = '$l_s$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if k == indxfixpbgalsour:
                    namefixp[k] = 'bgalsour'
                    strgfixp[k] = '$b_s$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if k == indxfixpspecsour:
                    namefixp[k] = 'specsour'
                    strgfixp[k] = '$f_s$'
                    scalfixp[k] = 'logt'
                if k == indxfixpsizesour:
                    namefixp[k] = 'sizesour'
                    strgfixp[k] = '$a_s$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                if k == indxfixpellpsour:
                    namefixp[k] = 'ellpsour'
                    strgfixp[k] = r'$\epsilon_s$'
                    scalfixp[k] = 'self'
                if k == indxfixpanglsour:
                    namefixp[k] = 'anglsour'
                    strgfixp[k] = r'$\phi_s$'
                    scalfixp[k] = 'self'
                if k == indxfixplgalhost:
                    namefixp[k] = 'lgalhost'
                    strgfixp[k] = '$l_h$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if k == indxfixpbgalhost:
                    namefixp[k] = 'bgalhost'
                    strgfixp[k] = '$b_h$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if k == indxfixpspechost:
                    namefixp[k] = 'spechost'
                    strgfixp[k] = '$f_h$'
                    scalfixp[k] = 'logt'
                if k == indxfixpsizehost:
                    namefixp[k] = 'sizehost'
                    strgfixp[k] = '$a_h$'
                    scalfixp[k] = 'logt'
                if k == indxfixpbeinhost:
                    namefixp[k] = 'beinhost'
                    strgfixp[k] = r'$\theta_{E,h}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                if k == indxfixpellphost:
                    namefixp[k] = 'ellphost'
                    strgfixp[k] = r'$\epsilon_h$'
                    scalfixp[k] = 'self'
                if k == indxfixpanglhost:
                    namefixp[k] = 'anglhost'
                    strgfixp[k] = r'$\phi_h$'
                    scalfixp[k] = 'self'
                if k == indxfixpsherhost:
                    namefixp[k] = 'sherhost'
                    strgfixp[k] = r'$\gamma_e$'
                    scalfixp[k] = 'self'
                if k == indxfixpsanghost:
                    namefixp[k] = 'sanghost'
                    strgfixp[k] = r'$\phi_{\gamma}$'
                    scalfixp[k] = 'self'

        if scalfixp[k] == 'pois' or scalfixp[k] == 'self' or scalfixp[k] == 'logt' or scalfixp[k] == 'atan':
            
            if namefixp[k][:-1].endswith('pop'):
                l = int(namefixp[k][-1])
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k][:-4])[l]
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k][:-4])[l]
            elif namefixp[k][-1].isdigit():
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k][:-1])
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k][:-1])
            else:
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k])
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k])
        
        if scalfixp[k] == 'gaus' or scalfixp[k] == 'eerr':
            if gdat.psfninfoprio:
                meanfixp[k] = getattr(gdat, 'meanpsfp')[k-indxfixppsfp[0]]
                stdvfixp[k] = getattr(gdat, 'stdvpsfp')[k-indxfixppsfp[0]]
            else:
                if namefixp[k][:-1].endswith('pop'):
                    l = int(namefixp[k][-1])
                    meanfixp[k] = getattr(gdat, 'mean' + namefixp[k][:-4])[l]
                    stdvfixp[k] = getattr(gdat, 'stdv' + namefixp[k][:-4])[l]
                else:
                    meanfixp[k] = getattr(gdat, 'mean' + namefixp[k])
                    stdvfixp[k] = getattr(gdat, 'stdv' + namefixp[k])
        if scalfixp[k] == 'self':
            factfixp[k] = maxmfixp[k] - minmfixp[k]
        if scalfixp[k] == 'logt':
            factfixp[k] = log(maxmfixp[k] / minmfixp[k])
        if scalfixp[k] == 'atan':
            factfixp[k] = arctan(maxmfixp[k]) - arctan(minmfixp[k])
        if scalfixp[k] == 'gaus':
            minmfixp[k] = meanfixp[k] - 3. * stdvfixp[k]
            maxmfixp[k] = meanfixp[k] + 3. * stdvfixp[k]
        if scalfixp[k] == 'eerr':
            cdfnminmfixp[k], cdfndifffixp[k] = retr_eerrnorm(minmfixp[k], maxmfixp[k], meanfixp[k], stdvfixp[k])
        
        if k in indxfixpfluxdistbrek:
            strgfixpunit[k] = strgfixp[k] + ' [%s]' % gdat.strgfluxunit
        elif k in indxfixpsigc or k in indxfixpsigt:
            strgfixpunit[k] = strgfixp[k] + ' [%s]' % gdat.strganglunit
        else:
            strgfixpunit[k] = strgfixp[k]

    for attr, valu in locals().iteritems():
        if attr != 'gdat' and '__' not in attr and not attr.endswith('temp') and attr != 'cntr':
            setattr(gdat, strgpara + attr, valu)


def defn_defa(gdat, valu, strg, strgpara=''):

    numbpopl = getattr(gdat, strgpara + 'numbpopl')
    varb = getattr(gdat, strgpara + strg)
    if varb == None:
        valutemp = zeros(numbpopl) + valu
    else:
        valutemp = varb
    
    temp = getattr(gdat, strgpara + 'fixp')
    temp[getattr(gdat, strgpara + 'indxfixp' + strg)] = valutemp
    setattr(gdat, strgpara + 'fixp', temp)


def setp_varbfull(gdat, strgpara, listfeat, typelimt='minmmaxm', numbpopl=None):
    
    numbfeat = len(listfeat)

    if numbpopl != None:
        listfeattemp = []
        listfeattemp.append(zeros(numbpopl) + listfeat[0])
        listfeattemp.append(zeros(numbpopl) + listfeat[1])
    else:
        listfeattemp = listfeat

    if typelimt == 'minmmaxm':
        minmpara = getattr(gdat, 'minm' + strgpara)
        maxmpara = getattr(gdat, 'minm' + strgpara)
        if minmpara == None:
            setattr(gdat, 'minm' + strgpara, listfeattemp[0])
        if maxmpara == None:
            setattr(gdat, 'maxm' + strgpara, listfeattemp[1])
    else:
        meanpara = getattr(gdat, 'mean' + strgpara)
        stdvpara = getattr(gdat, 'stdv' + strgpara)
        if meanmpara == None:
            setattr(gdat, 'mean' + strgpara, listfeattemp[0])
        if stdvpara == None:
            setattr(gdat, 'stdv' + strgpara, listfeattemp[1])
 

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
    indxoaxipnts = digitize(oaxi[0], gdat.binsoaxiopen) - 1
   
    return indxoaxipnts


def init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=None, indxevttplot=None, indxpoplplot=None):

    if strg == 'this' or gdatmodi != None:
        pathfold = gdat.pathfram
    elif strg == 'true' or strg == '':
        pathfold = gdat.pathinit
    elif strg == 'post':
        pathfold = gdat.pathpost
    
    figr, axis = plt.subplots(figsize=(gdat.sizeimag, gdat.sizeimag))
    
    if indxenerplot == None:
        strgener = 'A'
    else:
        strgener = '%d' % gdat.indxenerincl[indxenerplot]
    
    if indxevttplot == None:
        strgevtt = 'A'
    else:
        strgevtt = '%d' % gdat.indxevttincl[indxevttplot]
    
    if indxpoplplot == None:
        strgpopl = 'A'
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


def retr_imag(gdat, axis, maps, thisindxener=None, thisindxevtt=None, cmap='Reds', mean=False, vmin=None, vmax=None, scal=None):
    
    if scal == None:
        scal = gdat.scalmaps

    if vmin == None and vmax != None:
        vmin = -vmax
    
    draw_frambndr(gdat, axis)
   
    # filter the map
    if thisindxevtt == None:
        if thisindxener != None:
            if mean:
                maps = sum(maps[thisindxener, :, :] * gdat.expo[thisindxener, :, :], axis=1) / sum(gdat.expo[thisindxener, :, :], axis=1)
            else:
                maps = sum(maps[thisindxener, :, :], axis=1)
    else:
        maps = maps[thisindxener, :, thisindxevtt]
    
    # project the map to 2D
    if gdat.pixltype == 'heal':
        maps = tdpy.util.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                            minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                            minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
    if gdat.pixltype == 'cart':
        mapstemp = empty(gdat.numbsidecart**2)
        mapstemp[gdat.indxpixlrofi] = maps
        maps = mapstemp.reshape((gdat.numbsidecart, gdat.numbsidecart)).T
   
    if gdat.numbener > 1 and thisindxevtt == None and thisindxener == None:
        # plot the color of the map
        mapstemp = sum(maps, 2)
        mapstemp = maps[0, :] / maps[-1, :]
        mapstemp /= amax(mapstemp)
        mapsoutp = zeros((gdat.numbpixl, 3))
        mapsoutp[0, :] = mapstemp
        mapsoutp[2, :] = 1. - mapstemp
        maps = mapsoutp
    else:
        # rescale the map
        if scal == 'asnh':
            maps = arcsinh(maps)
    
    imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest', vmin=vmin, vmax=vmax, alpha=gdat.alphmaps)
    
    return imag


def make_cbar(gdat, axis, imag, indxenerplot=None, tick=None, labl=None):

    # make a color bar
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05, aspect=15)
    if tick != None and labl != None:
        cbar.set_ticks(tick)
        cbar.set_ticklabels(labl)
    
    return cbar


def make_catllabl(gdat, axis):

    axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, label='Model PS', marker='+', linewidth=2, color='b')
    
    if gdat.trueinfo:
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablhits, marker='x', linewidth=2, color='g')
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablbias, marker='*', linewidth=2, color='g', facecolor='none')
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, facecolor='none', \
                                                                                                label=gdat.truelablmiss, marker='o', linewidth=2, color='g')
    if gdat.pntstype == 'lens':
        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='Model Source', marker='<', linewidth=2, color='b')

        axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='Model Host', marker='s', linewidth=2, color='b')
        if gdat.trueinfo:
            axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='%s Source' % gdat.truelabl, marker='>', linewidth=2, color='g')
        
            axis.scatter(gdat.anglfact * gdat.maxmgang * 5., gdat.anglfact * gdat.maxmgang * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='%s Host' % gdat.truelabl, marker='D', linewidth=2, color='g')
        
    axis.legend(bbox_to_anchor=[0.5, 1.1], loc='center', ncol=4)
        

def supr_fram(gdat, gdatmodi, axis, indxpoplplot=None, trueonly=False):

    # true catalog
    if gdat.trueinfo:
       
        if indxpoplplot == None:
            indxpoplplot = gdat.indxpopl
        else:
            indxpoplplot = [indxpoplplot]

        for l in indxpoplplot:
            ## get the true catalog
            if gdat.numbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist, :].flatten())
                lgal = copy(gdat.truelgal[l])
                bgal = copy(gdat.truebgal[l])
                numbpnts = int(gdat.truefixp[gdat.indxfixpnumbpnts][l])
                
                if gdatmodi != None and not trueonly:
                    
                    ## associations
                    ### missed
                    indx = gdatmodi.trueindxpntsassc[l].miss
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablmiss, facecolor='none', \
                                                                                                                                marker='o', linewidth=2, color='g')
                    
                    ### biased
                    indx = gdatmodi.trueindxpntsassc[l].bias[gdat.indxenerfluxdist]
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablbias, marker='*', linewidth=2, color='g', facecolor='none')
                    
                    ### hit
                    indx = gdatmodi.trueindxpntsassc[l].hits[gdat.indxenerfluxdist]
                    
                else:
                    indx = arange(lgal.size)
                
                axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                        label=gdat.truelablhits, marker='x', linewidth=2, color='g')
            
            if gdat.pntstype == 'lens':
               
                ## host
                axis.scatter(gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalhost], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalhost], \
                                                                            alpha=gdat.alphpnts, label=gdat.truelablhits, s=300, marker='D', linewidth=2, color='g')
                axis.add_patch(plt.Circle((gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalhost], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalhost]), \
                                                                                gdat.fluxfactplot * gdat.truefixp[gdat.trueindxfixpbeinhost], edgecolor='g', facecolor='none', lw=2))
                
                ## source
                axis.scatter(gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalsour], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalsour], \
                                                                                    alpha=gdat.alphpnts, label=gdat.truelablhits, s=300, marker='>', linewidth=2, color='g')
    
                if gdat.numbtrap > 0:
                    ## subhalos
                    for k in range(lgal.size):
                        axis.add_patch(plt.Circle((gdat.anglfact * lgal[k], gdat.anglfact * bgal[k]), \
                                        gdat.fluxfactplot * gdat.truespec[l][0, gdat.indxenerfluxdist[0], k], edgecolor='g', facecolor='none', lw=2))
    
            ## annotate
            if gdat.anotcatl:
                for a in range(numbpnts):
                    strg = ''
                    if gdat.truestrg[l][a] != None:
                        strg += '%s ' % gdat.truestrg[l][a]
                    if gdat.truestrgassc[l][a] != None:
                        strg += '%s ' % gdat.truestrgassc[l][a]
                    if gdat.truestrgclss[l][a] != None:
                        strg += '%s ' % gdat.truestrgclss[l][a]
                    if strg != '':
                        axis.text(gdat.anglfact * gdat.truelgal[l][a], gdat.anglfact * gdat.truebgal[l][a] - gdat.offstext, strg, \
                                                                                                                            ha='center', va='center', color='g', fontsize=6)
    
        # model catalog
        if gdatmodi != None and not trueonly:
            if gdat.numbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]])
                lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l]]
                bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l]]
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphpnts, label='Sample', marker='+', linewidth=2, color='b')
            if gdat.pntstype == 'lens':
                
                ## source
                lgalsour = gdatmodi.thissampvarb[gdat.indxfixplgalsour]
                bgalsour = gdatmodi.thissampvarb[gdat.indxfixpbgalsour]
                axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, alpha=gdat.alphpnts, label='Model Source', s=300, marker='<', linewidth=2, color='b')
    
                ## host
                lgalhost = gdatmodi.thissampvarb[gdat.indxfixplgalhost]
                bgalhost = gdatmodi.thissampvarb[gdat.indxfixpbgalhost]
                beinhost = gdatmodi.thissampvarb[gdat.indxfixpbeinhost]
                axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, alpha=gdat.alphpnts, label='Model Host', s=300, marker='s', linewidth=2, color='b')
                axis.add_patch(plt.Circle((gdat.anglfact * lgalhost, gdat.anglfact * bgalhost), gdat.fluxfactplot * beinhost, edgecolor='b', facecolor='none', lw=2, ls='--'))
                
                # subhalos
                if gdat.numbtrap > 0:
                    for k in range(lgal.size):
                        axis.add_artist(plt.Circle((gdat.anglfact * lgal[k], gdat.anglfact * bgal[k]), \
                                    gdat.fluxfactplot * gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][0, k]], edgecolor='b', facecolor='none', ls='--', lw=2))
    
    
def retr_indxpsfp(gdat, psfntype, varioaxi):

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

    numbpsfptotl = numbpsfpform + numbpsfpoaxi
    
    if varioaxi:
        indxpsfpoaxinorm = numbpsfpform + numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)
        indxpsfpoaxiindx = numbpsfpform + numbpsfptotl * arange(gdat.numbener * gdat.numbevtt) + 1
    else:
        indxpsfpoaxinorm = []
        indxpsfpoaxiindx = []

    return numbpsfpform, numbpsfpoaxi, numbpsfptotl, indxpsfpoaxinorm, indxpsfpoaxiindx


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
    
    if gdat.pixltype == 'heal':
        angldist = hp.rotator.angdist(dir1, dir2)
    else:
        angldist = sqrt((dir1[0, :] - dir2[0])**2 + (dir1[1, :] - dir2[1])**2)

    return angldist


def corr_catl(gdat, gdatmodi, thisindxpopl, modllgal, modlbgal, modlspec, modldeflpnts=None, metrtype='dist'):

    trueindxpntsassc = tdpy.util.gdatstrt()
    trueindxpntsassc.miss = []
    trueindxpntsassc.bias = [[] for i in gdat.indxener]
    trueindxpntsassc.hits = [[] for i in gdat.indxener]
    trueindxpntsassc.mult = []
        
    indxmodlpnts = zeros_like(gdat.truelgal[thisindxpopl], dtype=int) - 1
    specassc = zeros((gdat.numbener, gdat.truefixp[gdat.indxfixpnumbpnts[thisindxpopl]]), dtype=float)
    numbassc = zeros_like(gdat.truelgal[thisindxpopl], dtype=int)
    metrassc = zeros_like(gdat.truelgal[thisindxpopl]) + 3 * gdat.maxmgang
    
    if metrtype == 'dist':
        dir1 = array([gdat.truelgal[thisindxpopl], gdat.truebgal[thisindxpopl]])

    for k in range(modllgal.size):
       
        # determine which true PSs satisfy the match criterion
        if metrtype == 'dist':
            dir2 = array([modllgal[k], modlbgal[k]])
            metr = retr_angldist(gdat, dir1, dir2)
            trueindxpntstemp = where(metr < gdat.anglassc)[0]
        if metrtype == 'defl':
            metr = sum(sum(sum(gdat.truedeflpnts * modldeflpnts[:, :, :, k], 0, 0, 0))) / gdat.numbpixl / 2. / gdat.pixlsize**2
            trueindxpntstemp = where(metr > 1.)[0]
        
        if trueindxpntstemp.size > 0:
            
            # if there are multiple associated true PS, sort them
            indx = argsort(metr[trueindxpntstemp])
            metr = metr[trueindxpntstemp][indx]
            trueindxpntstemp = trueindxpntstemp[indx]
                
            # store the index of the model PS
            numbassc[trueindxpntstemp[0]] += 1
            if metr[0] < metrassc[trueindxpntstemp[0]]:
                specassc[:, trueindxpntstemp[0]] = modlspec[:, k]
                metrassc[trueindxpntstemp[0]] = metr[0]
                indxmodlpnts[trueindxpntstemp[0]] = k

    # get the flux limit that delineates the biased associations and hits 
    fluxbias = empty((2, gdat.numbener, gdat.truefixp[gdat.indxfixpnumbpnts[thisindxpopl]]))
    for i in gdat.indxener:
        fluxbias[:, i, :] = retr_fluxbias(gdat, gdat.truespec[thisindxpopl][0, i, :], i)

    # divide associations into subgroups
    for k in range(gdat.truefixp[gdat.indxfixpnumbpnts[thisindxpopl]].astype(int)):
        if numbassc[k] == 0:
            trueindxpntsassc.miss.append(k)
        else:
            if numbassc[k] > 1:
                trueindxpntsassc.mult.append(k)
    
            ## check whether the flux of the associated model point source matches well with the flux of the deterministic point source
            for i in gdat.indxener:
                boolbias = specassc[i, k] > fluxbias[1, i, k] or specassc[i, k] < fluxbias[0, i, k]
                if boolbias:
                    trueindxpntsassc.bias[i].append(k)
                else:
                    trueindxpntsassc.hits[i].append(k)
   
    if False and gdat.verbtype > 1:
        print 'Correlating catalogs...'
        print 'thisindxpopl'
        print thisindxpopl
        print 'trueindxpntsassc.hits'
        print trueindxpntsassc.hits
        print 'trueindxpntsassc.bias'
        print trueindxpntsassc.bias
        print 'trueindxpntsassc.mult'
        print trueindxpntsassc.mult
        print 'trueindxpntsassc.miss'
        print trueindxpntsassc.miss
        print 

    return indxmodlpnts, trueindxpntsassc


def pert_llik(gdat, gdatmodi, indxparapert, stdvparapert):

    numbpert = indxparapert.size 
    gdatmodi.drmcsamp[:, 1] = gdatmodi.drmcsamp[:, 0]
    for k in range(numbpert):
        gdatmodi.drmcsamp[indxparapert[k], 1] += stdvparapert[k]
    
    gdatmodi.nextsampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.drmcsamp[:, 1], 'next')
    
    lpospert = retr_negalpos(gdat, gdatmodi)
    
    return lpospert


def retr_defl(gdat, lgal, bgal, bein, ellp, angl, sher, sang, rcor):
        
    # translate the grid
    lgaltran = gdat.lgalgridcart - lgal
    bgaltran = gdat.bgalgridcart - bgal
    
    # rotate the grid
    lgalrttr = cos(angl) * lgaltran - sin(angl) * bgaltran
    bgalrttr = sin(angl) * lgaltran + cos(angl) * bgaltran
    
    if ellp > 1e-4:
        factflat = (1. - ellp)**2
        factrcor = sqrt(factflat * (rcor**2 + lgalrttr**2) + bgalrttr**2)
        facteccc = sqrt(1. - factflat)
        factbein = (bein * (1. - ellp) / facteccc)
        defllgalrttr = factbein *  arctan(facteccc * lgalrttr / (factrcor + rcor))
        deflbgalrttr = factbein * arctanh(facteccc * bgalrttr / (factrcor + factflat * rcor))
    else:
        rint = sqrt(lgalrttr**2 + bgalrttr**2 + rcor**2)
        
        defllgalrttr = bein * lgalrttr / (rint + rcor) 
        deflbgalrttr = bein * bgalrttr / (rint + rcor)
    
    # totate back vector to original basis
    defllgal =  cos(angl) * defllgalrttr + sin(angl) * deflbgalrttr
    deflbgal = -sin(angl) * defllgalrttr + cos(angl) * deflbgalrttr
    
    # external shear
    factcosi = sher * cos(2. * sang)
    factsine = sher * cos(2. * sang)
    defllgal += factcosi * gdat.lgalgridcart + factsine * gdat.bgalgridcart
    deflbgal += factsine * gdat.lgalgridcart - factcosi * gdat.bgalgridcart
    
    return dstack((defllgal, deflbgal))
   

def retr_negalpos(gdat, gdatmodi):
   
    proc_samp(gdat, gdatmodi, 'next')
    
    return -gdatmodi.nextlliktotl - gdatmodi.nextlpritotl


def retr_lpridist(gdat, gdatmodi, flux, sind, curv, expo, lpri, sampvarb):
    
    for l in gdat.indxpopl:
        if gdat.fluxdisttype == 'powr':
            fluxdistslop = sampvarb[gdat.indxfixpfluxdistslop[l]]
            lpriflux = sum(log(pdfn_flux_powr(gdat, flux, fluxdistslop)))
        if gdat.fluxdisttype == 'brok':
            fluxdistbrek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
            fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
            fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
            lpriflux = sum(log(pdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)))
        lpri[1+gdat.numbpopl+l] = lpriflux

        if gdat.numbener > 1:
            sinddistmean = sampvarb[gdat.indxfixpsinddistmean[l]]
            sinddiststdv = sampvarb[gdat.indxfixpsinddiststdv[l]]
            lpri[1+2*gdat.numbpopl+l] = sum(lpdf_gaus(sind, gdatmodi.tempsinddistmean[l], gdatmodi.tempsinddiststdv[l])) 
            if gdat.spectype[l] == 'curv':
                curvdistmean[l] = sampvarb[gdat.indxfixpcurvdistmean[l]]
                curvdiststdv[l] = sampvarb[gdat.indxfixpcurvdiststdv[l]]
                lpri[1+3*gdat.numbpopl+l] = sum(lpdf_gaus(sind, curvdistmean[l], curvdiststdv[l])) 
            if gdat.spectype[l] == 'expo':
                expodistmean[l] = sampvarb[gdat.indxfixpexpodistmean[l]]
                expodiststdv[l] = sampvarb[gdat.indxfixpexpodiststdv[l]]
                lpri[1+3*gdat.numbpopl+l] = sum(lpdf_gaus(sind, expodistmean[l], expodiststdv[l])) 


def proc_samp(gdat, gdatmodi, strg, raww=False):

    if gdat.verbtype > 1:
        print 'proc_samp()'
        print 'strg'
        print strg
        print

    if gdatmodi != None:
        gdatobjt = gdatmodi
    else:
        gdatobjt = gdat

    if strg == 'true':
        strgtype = 'true'
        spectype = gdat.truespectype
    else:
        strgtype = ''
        spectype = gdat.spectype
    
    # grab the sample vector
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
        sizehost = sampvarb[getattr(gdat, strgtype + 'indxfixpsizehost')]
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

    indxpntsfull = list(getattr(gdatobjt, strg + 'indxpntsfull'))
    indxsamplgal, indxsampbgal, indxsampflux, indxsampspec, indxsampsind, indxsampcurv, indxsampexpo, indxsampcompcolr = retr_indx(gdat, indxpntsfull, spectype)
    numbpnts = getattr(gdatobjt, strg + 'sampvarb')[gdat.indxfixpnumbpnts].astype(int)
    lgal = []
    bgal = []
    flux = []
    spec = []
    numbpopl = numbpnts.size
    for l in range(numbpopl):
        lgal.append(sampvarb[indxsamplgal[l]])
        bgal.append(sampvarb[indxsampbgal[l]])
        flux.append(sampvarb[indxsampflux[l]])
        spec.append(sampvarb[indxsampspec[l]])
    lgalconc = concatenate(lgal)
    bgalconc = concatenate(bgal)
    specconc = concatenate(spec, axis=1)
    if gdat.numbener > 1:
        indxsampspec = getattr(gdatobjt, strg + 'indxsampspec')
        sind = [[] for l in gdat.indxpopl]
        curv = [[] for l in gdat.indxpopl]
        expo = [[] for l in gdat.indxpopl]
        for l in range(numbpopl):
            sind[l] = sampvarb[indxsampsind[l]]
            if gdat.spectype == 'curv':
                curv[l] = sampvarb[indxsampcurv[l]]
            if gdat.spectype == 'expo':
                expo[l] = sampvarb[indxsampexpo[l]]
        
    if strg == 'next' and gdat.verbtype > 1:
        setattr(gdatobjt, strg + 'indxsamplgal', indxsamplgal)
        setattr(gdatobjt, strg + 'indxsampbgal', indxsampbgal)
        setattr(gdatobjt, strg + 'indxsampflux', indxsampflux)
        setattr(gdatobjt, strg + 'indxsampspec', indxsampspec)
        setattr(gdatobjt, strg + 'indxsampsind', indxsampsind)
        setattr(gdatobjt, strg + 'indxsampcurv', indxsampcurv)
        setattr(gdatobjt, strg + 'indxsampexpo', indxsampexpo)
        show_samp(gdat, gdatmodi)
    
    numbpntsconc = lgalconc.size
    
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
    
    if strg != 'true':
        setattr(gdatobjt, strg + 'lgal', lgal)
        setattr(gdatobjt, strg + 'bgal', bgal)
        setattr(gdatobjt, strg + 'flux', flux)
        setattr(gdatobjt, strg + 'spec', spec)
        if gdat.numbener > 1:
            setattr(gdatobjt, strg + 'sind', sind)
            setattr(gdatobjt, strg + 'curv', curv)
            setattr(gdatobjt, strg + 'expo', expo)
    
    bacp = sampvarb[getattr(gdat, 'indxfixpbacp')]
    
    if gdat.pntstype == 'lens':
        
        ## components
        if numbpntsconc > 0 and not raww:
            for k in range(numbpntsconc):
                # temp -- fix truncation
                defl += retr_defl(gdat, lgalconc[k], bgalconc[k], specconc[0, k], 0., 0., 0., 0., 0.)
                        
        # lensed image
        lensflux = retr_mapsraww(gdat, gdat.lgalgridcart - defl[:, :, 0], gdat.bgalgridcart - defl[:, :, 1], bacp, lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
        
        # host emission
        hostfluxmaps = retr_mapssers(gdat, gdat.lgalgridcart, gdat.bgalgridcart, lgalhost, bgalhost, spechost, sizehost, ellphost, anglhost)
        
        # total emission
        modlfluxuncv = lensflux + bacp * gdat.backfluxcart[0] + hostfluxmaps
        
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
        
        if False:
            print 'defl'
            summgene(defl)
            print 'lensflux'
            summgene(lensflux)
            print 'hostfluxmaps'
            summgene(hostfluxmaps)
            print 'modlflux'
            summgene(modlflux)
            print

    if gdat.pntstype == 'lght':
        
        ### PS flux map
        pntsflux = retr_pntsflux(gdat, lgalconc, bgalconc, specconc, psfnintp, varioaxi, evalcirc=gdat.evalcirc)
        setattr(gdatobjt, strg + 'pntsflux', pntsflux)
        
        ### model flux map
        modlflux = retr_mapslght(gdat, bacp, pntsflux, gdat.indxcube)
   
    ### count map
    if gdat.pixltype != 'unbd':
        modlcnts = retr_cntsmaps(gdat, modlflux)
        setattr(gdatobjt, strg + 'modlcnts', modlcnts)
        if gdat.verbtype > 1:
            print 'modlcnts'
            summgene(modlcnts)
    
    meanpnts = sampvarb[gdat.indxfixpmeanpnts]
    if strg != 'true':
        
        ### log-prior
        if gdat.numbtrap > 0:
            for l in gdat.indxpopl:
                lpri = zeros(gdat.numblpri)
                lpri[0] = -0.5 * gdat.priofactdoff * gdat.numbcompcolr[l] * numbpnts[l]
                lpri[1+l] = retr_probpois(numbpnts[l], meanpnts[l])
                retr_lpridist(gdat, gdatmodi, flux[l], sind[l], curv[l], expo[l], lpri, sampvarb)
            lpritotl = sum(lpri)
            
            if strg == 'next':
                prevlpri = copy(lpri)
                retr_lpridist(gdat, gdatmodi, flux[l], sind[l], curv[l], expo[l], prevlpri, gdatmodi.thissampvarb)
                deltlpri = lpritotl - sum(prevlpri)
                setattr(gdatmodi, 'deltlpri', deltlpri) 
            
        ### log-likelihood
        if gdat.pixltype == 'unbd':
            gdatmodi.templlik = gdat.numbdatasamp * log(modlfluxtotl) - modlfluxtotl + log(modlflux)
        else:
            if gdat.liketype == 'pois':
                gdatmodi.templlik = gdat.datacnts * log(modlcnts) - modlcnts
            if gdat.liketype == 'gaus':
                gdatmodi.templlik = -0.5 * (gdat.datacnts - modlcnts)**2 / gdat.datacnts
           
        gdatmodi.templliktotl = sum(gdatmodi.templlik)
        if strg == 'next':
            deltllik = gdatmodi.templliktotl - gdatmodi.thislliktotl
            setattr(gdatmodi, 'deltllik', deltllik)
            # temp
            laccfact = 0.
            setattr(gdatmodi, 'laccfact', laccfact)
       
        setattr(gdatobjt, strg + 'lpritotl', lpritotl) 
        setattr(gdatobjt, strg + 'lliktotl', gdatmodi.templliktotl) 
        
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
        
        if not raww:
            retr_datatick(gdat)
        else:
            return

    ## tertiary variables that are not needed for likelihood evaluation
    if strg != 'true':
        setattr(gdatmodi, strg + 'lpri', lpri)
        setattr(gdatmodi, strg + 'llik', gdatmodi.templlik) 
    if strg != 'next':
        if strg == 'this':
            lpostotl = lpritotl + gdatmodi.templliktotl
            setattr(gdatobjt, strg + 'lpostotl', lpostotl) 
       
        setattr(gdatobjt, strg + 'psfp', psfp)
        if gdat.pixltype != 'unbd':
            resicnts = gdat.datacnts - modlcnts
            setattr(gdatobjt, strg + 'resicnts', resicnts)
            
        ## prior on the flux distribution
        if gdat.numbtrap > 0:
            fluxhistprio = empty((numbpopl, gdat.numbfluxplot))
            for l in range(numbpopl):
                if gdat.fluxdisttype == 'powr':
                    fluxdistslop = sampvarb[gdat.indxfixpfluxdistslop[l]]
                    fluxhistprio[l, :] = meanpnts[l] * pdfn_flux_powr(gdat, gdat.meanfluxplot, fluxdistslop) * gdat.deltfluxplot
                if gdat.fluxdisttype == 'brok':
                    fluxdistbrek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
                    fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
                    fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
                    fluxhistprio[l, :] = meanpnts[l] * pdfn_flux_brok(gdat, gdat.meanfluxplot, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr) * gdat.deltfluxplot
            setattr(gdatobjt, strg + 'fluxhistprio', fluxhistprio)
        
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
                pntscnts = retr_cntsmaps(gdat, pntsflux)
                errrcnts = pntscnts - temppntscnts
                indxcubegood = where(temppntscnts > 1e-10)
                setattr(gdatobjt, strg + 'errrcnts', errrcnts)
                if False and amax(fabs(errrcnts)) > 0.1:
                    raise Exception('Approximation error in calculating the PS flux map is above the tolerance level.')

        if gdat.pntstype == 'lens':

            
            lenscnts = retr_cntsmaps(gdat, lensflux, cart=True)
            hostcntsmaps = retr_cntsmaps(gdat, hostfluxmaps, cart=True)
            
            setattr(gdatobjt, strg + 'psfnkern', psfnkern)
            setattr(gdatobjt, strg + 'modlfluxuncv', modlfluxuncv)
            setattr(gdatobjt, strg + 'lenscnts', lenscnts)
            setattr(gdatobjt, strg + 'hostcntsmaps', hostcntsmaps)

            ### deflection
            #### number of deflection components
            numbdeflpnts = min(3, numbpntsconc)
            if numbpntsconc > 0:
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
                gdatmodi.thisdeflcomp = 1. - sum(gdatmodi.thisdefl * gdat.truedefl, 2) / sqrt(gdatmodi.thisdefl[:, :, 0]**2 + gdatmodi.thisdefl[:, :, 1]**2) / \
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
                if gdat.fluxdisttype == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[l]]
                    gdatmodi.thislprinorm -= log(1. + fluxdistslop**2)
                    gdatmodi.thislprinorm += sum(log(pdfn_flux_powr(gdat, flux, fluxdistslop)))
                if gdat.fluxdisttype == 'brok':
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


def retr_mapsgaus(gdat, lgal, bgal, spec, size, ellp, angl):
    
    rttrmatr = array([[cos(angl), -sin(angl)], [sin(angl), cos(angl)]])
    icovmatr = array([[1. / ((1. - ellp) * size)**2, 0.], [0., 1. / size**2]])

    posi = array([lgalgrid - lgal, bgalgrid - bgal])
    mapsgaus = flux * exp(-0.5 * sum(posi * tensordot(self.icovmatr, posi, (1,0)), 0)) / size**2 / (1. - ellp)
        
    return mapsgaus


def retr_mapssers(gdat, lgalgrid, bgalgrid, lgal, bgal, spec, size, ellp, angl):
    
    lgalrttr = (1. - ellp) * (cos(angl) * (lgalgrid - lgal) - sin(angl) * (bgalgrid - bgal))
    bgalrttr = sin(angl) * (lgalgrid - lgal) + cos(angl) * (bgalgrid - bgal) 

    # temp
    #mapssers = spec[:, None, None] * exp(-7.67 * ((sqrt(lgalrttr[None, :, :]**2 + bgalrttr[None, :, :]**2) / size)**0.25 - 1.))

    mapssers = 1e-3 * spec * exp(-7.67 * ((sqrt(lgalrttr**2 + bgalrttr**2) / size)**0.25 - 1.)) / gdat.apix
    mapssers = mapssers[None, :, :, None]

    return mapssers

    
def retr_mapsraww(gdat, lgalgrid, bgalgrid, bacp, lgal, bgal, spec, size, ellp, angl):
    
    # source emission
    mapssour = retr_mapssers(gdat, lgalgrid, bgalgrid, lgal, bgal, spec, size, ellp, angl)
    
    # total unlensed emission
    mapsraww = mapssour
        
    return mapsraww


