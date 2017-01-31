# common imports
from __init__ import *

class interp1d_pick:
    """ class wrapper for piecewise linear function
    """
    def __init__(self, xi, yi, **kwargs):
        self.xi = xi
        self.yi = yi
        self.args = kwargs
        self.f = interp1d(xi, yi, **kwargs)

    def __call__(self, xnew):
        return self.f(xnew)

    def __getstate__(self):
        return self.xi, self.yi, self.args

    def __setstate__(self, state):
        self.f = interp1d(state[0], state[1], **state[2])


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
                    wdthtemp = interp1d_pick(psfntemp[indxanglgood], gdat.binsangl[indxanglgood])(intpwdth)
                if varioaxi:
                    wdth[i, m, p] = wdthtemp
                else:
                    wdth[i, m] = wdthtemp
                        
    return wdth


def retr_spec(gdat, flux, sind, curv, expo, spectype, plot=False):
    
    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if plot:
            meanener = gdat.meanenerplot
        else:
            meanener = gdat.meanener

        if spectype == 'powr':
            spec = flux[None, :] * (meanener / gdat.enerfluxdist)[:, None]**(-sind[None, :])
        if spectype == 'curv':
            spec = flux[None, :] * meanener[:, None]**(-sind[None, :] - gdat.factlogtenerpivt[:, None] * curv[None, :])
        if spectype == 'expo':
            spec = flux[None, :] * (meanener / gdat.enerfluxdist)[:, None]**(-sind[None, :]) * exp(-(meanener - gdat.enerfluxdist)[:, None] / expo[None, :])
    
    return spec


def retr_indx(gdat, indxpntsfull, spectype):
    
    indxsamplgal = [[] for l in gdat.indxpopl]
    indxsampbgal = [[] for l in gdat.indxpopl]
    indxsampflux = [[] for l in gdat.indxpopl]
    indxsampsind = [[] for l in gdat.indxpopl]
    indxsampcurv = [[] for l in gdat.indxpopl]
    indxsampexpo = [[] for l in gdat.indxpopl]
    indxsampcomp = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            indxsamptemp = gdat.indxsampcomp[0] + gdat.numbtrapcuml[l] + array(indxpntsfull[l], dtype=int) * gdat.numbcomp[l]
            indxsamplgal[l] = indxsamptemp
            indxsampbgal[l] = indxsamptemp + 1
            indxsampflux[l] = indxsamptemp + 2
            if gdat.numbener > 1:
                indxsampsind[l] = indxsamptemp + 3
                if spectype[l] == 'curv':
                    indxsampcurv[l] = indxsamptemp + 4
                if spectype[l] == 'expo':
                    indxsampexpo[l] = indxsamptemp + 4
            indxsampcomp.append(repeat(indxsamptemp, gdat.numbcomp[l]) + tile(gdat.indxcomp[l], len(indxpntsfull[l])))
             
    return indxsamplgal, indxsampbgal, indxsampflux, indxsampsind, indxsampcurv, indxsampexpo, indxsampcomp


def retr_plotpath(gdat, gdatmodi, strg, strgplot):
    
    if gdatmodi == None:
        if strg == 'true':
            path = gdat.pathinit + strgplot + '.pdf'
        else:
            path = gdat.pathpost + strgplot + '.pdf'
    else:
        path = gdat.pathfram + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


def retr_pntsflux(gdat, lgal, bgal, spec, psfnintp, varioaxi, evalcirc):
    
    if gdat.verbtype > 1:
        print 'retr_pntsflux'

    numbpnts = lgal.size
    if gdat.pixltype == 'unbd':
        pntsfluxsing = zeros((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt, 2))
    else:
        pntsfluxsing = zeros((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    if False:
        print 'lgal'
        print lgal
        print 'bgal'
        print bgal
        print 'spec'
        print spec
    
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


def pdfn_brok(xdat, brek, sloplowr, slopuppr, minm, maxm):

    norm = 1. / (brek**sloplowr * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) / (1. - sloplowr) + \
                 brek**slopuppr * (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    pdfn = norm * (xdat / brek)**(-sloplowr)
    indx = where(xdat >= brek)[0]
    if indx.size > 0:
        pdfn[indx] = norm * (xdat[indx] / brek)**(-slopuppr)
        
    return pdfn


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


def cdfn_powr(flux, minm, maxm, slop):
        
    unit = (flux**(1. - slop) - minm**(1. - slop)) / (maxm**(1. - slop) - minm**(1. - slop))
        
    return unit


def icdf_powr(unit, minm, maxm, slop):

    para = (unit * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))

    return para


def icdf_expo(unit, maxm, scal):

    para = -scal * log(1. - unit * (1. - exp(-maxm / scal)))

    return para


def pdfn_expo(xdat, maxm, scal):

    if (xdat > maxm).any():
        pdfn = 0.
    else:
        pdfn = 1. / scal / (1. - exp(-maxm / scal)) * exp(-xdat / scal)

    return pdfn


def pdfn_dexp(xdat, maxm, scal):
    
    pdfn = 0.5 * pdfn_expo(fabs(xdat), maxm, scal)

    return pdfn


def pdfn_powr(xdat, slop, minm, maxm):
  
    norm = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop))
    
    pdfn = norm * xdat**(-slop)
    
    return pdfn


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


def pdfn_gaus(xdat, mean, stdv):
    
    pdfn = 1. / sqrt(2. * pi) / stdv * exp(-0.5 * ((xdat - mean) / stdv)**2)

    return pdfn


def cdfn_lgau(para, mean, stdv):
    
    cdfn = cdfn_gaus(log(para), log(mean), stdv)

    return cdfn


def icdf_lgau(paraunit, mean, stdv):
    
    icdf = exp(icdf_gaus(paraunit, log(mean), stdv))

    return icdf


def pdfn_lgau(xdat, mean, stdv):
    
    pdfn = pdfn_gaus(log(xdat), log(mean), stdv)

    return pdfn


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
            gdatmodi.thisindxprop = gdat.indxproptypebrth
        else:
            gdatmodi.thisindxprop = gdat.indxproptypedeth

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
            gdatmodi.thisindxprop = gdat.indxproptypesplt
        else:
            gdatmodi.thisindxprop = gdat.indxproptypemerg

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
                raise Exception('pixlcnvt went negative!')

    if gdat.pixltype == 'cart':
        indxlgcr = floor(gdat.numbsidecart * (lgal - gdat.minmlgaldata) / 2. / gdat.maxmgangdata).astype(int)
        indxbgcr = floor(gdat.numbsidecart * (bgal - gdat.minmbgaldata) / 2. / gdat.maxmgangdata).astype(int)

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
        if gdat.enerdiff:
            cntsmaps *= gdat.deltener[:, None, None, None]
    else:
        
        cntsmaps = fluxmaps * gdat.expo * gdat.apix
        if gdat.enerdiff:
            cntsmaps *= gdat.deltener[:, None, None]
        
    return cntsmaps


def retr_cntsbackfwhm(gdat, bacp, fwhm):

    varioaxi = len(fwhm.shape) == 3
    cntsbackfwhm = zeros_like(fwhm)
    for c in gdat.indxback:
        indxbacp = c * gdat.numbener + gdat.indxener
        if varioaxi:
            cntsback = bacp[indxbacp, None, None, None] * gdat.backflux[c][:, :, :, None] * gdat.expo[:, :, :, None] * \
                                                                                                gdat.deltener[:, None, None, None] * pi * fwhm[:, None, :, :]**2 / 4.
        else:
            cntsback = bacp[indxbacp, None, None] * gdat.backflux[c] * gdat.expo * pi * fwhm[:, None, :]**2 / 4.
            if gdat.enerdiff:
                cntsback *= gdat.deltener[:, None, None]
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
    
    # temp
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
    indxsamplgal, indxsampbgal, indxsampflux, indxsampsind, indxsampcurv, indxsampexpo, indxsampcomp = retr_indx(gdat, indxpntsfull, spectype) 
    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxfixpnumbpnts] = samp[gdat.indxfixpnumbpnts]
    
    for k in gdat.indxfixp:
        sampvarb[k] = icdf_fixp(gdat, '', samp[k], k)
    
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgang, 2. * gdat.maxmgang)
            sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgang, 2. * gdat.maxmgang) 
            
            if gdat.fluxdisttype[l] == 'powr':
                sampvarb[indxsampflux[l]] = icdf_flux_powr(samp[indxsampflux[l]], gdat.minmflux, gdat.maxmflux, sampvarb[gdat.indxfixpfluxdistslop[l]])
            if gdat.fluxdisttype[l] == 'brok':
                fluxunit = samp[indxsampflux[l]]
                fluxdistbrek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
                fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
                fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
                sampvarb[indxsampflux[l]] = icdf_flux_brok(fluxunit, gdat.minmflux, gdat.maxmflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
            
            if not isfinite(sampvarb[indxsampflux[l]]).all():
                print 'sampvarb[indxsampflux[l]]'
                print sampvarb[indxsampflux[l]]
                raise Exception('Flux is not finite.')
            
            if gdat.numbener > 1:
                sampvarb[indxsampsind[l]] = icdf_gaus(samp[indxsampsind[l]], sampvarb[gdat.indxfixpsinddistmean[l]], sampvarb[gdat.indxfixpsinddiststdv[l]])
                if gdat.spectype[l] == 'curv':
                    sampvarb[indxsampcurv[l]] = icdf_gaus(samp[indxsampcurv[l]], sampvarb[gdat.indxfixpcurvdistmean[l]], sampvarb[gdat.indxfixpcurvdiststdv[l]])
                if gdat.spectype[l] == 'expo':
                    sampvarb[indxsampexpo[l]] = icdf_gaus(samp[indxsampexpo[l]], sampvarb[gdat.indxfixpexpodistmean[l]], sampvarb[gdat.indxfixpexpodiststdv[l]])

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
            fermform[:, m, k] = interp1d_pick(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)
    
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
    gdatmodi.thissamp[gdatmodi.indxsampmodi] = gdatmodi.nextsamp[gdatmodi.indxsampmodi]
    
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

    gdat.exprsind = -log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :]) / log(gdat.meanener / gdat.enerfluxdist)
    
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
                                                                                            fgl3['Flux10000_100000']))[gdat.indxenerincl, :] / gdat.deltener[:, None]
    
    if False:
        print 'gdat.exprspec'
        for k in range(gdat.exprspec.shape[2]):
            print gdat.exprspec[:, 1, k].T
    
    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.exprspec[where(isfinite(gdat.exprspec) == False)] = 0.
   
    if False:
        print 'fgl3specstdvtemp0'
        print fgl3specstdvtemp[:, :, 0]
        print 'fgl3specstdvtemp1'
        print fgl3specstdvtemp[:, :, 1]

        print 'gdat.exprspec'
        for k in range(gdat.exprspec.shape[2]):
            print gdat.exprspec[:, 1, k].T

    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3strg = fgl3['Source_Name']
    fgl3strgclss = fgl3['CLASS1']
    fgl3strgassc = fgl3['ASSOC1']
    
    fgl3spectype = fgl3['SpectrumType']
    gdat.exprsind = fgl3['Spectral_Index']
    gdat.exprcurv = fgl3['beta']
    gdat.exprexpo = fgl3['Cutoff'] * 1e-3
    

def retr_rtag(gdat):
    
    rtag = '%d' % (gdat.numbswep)
        
    return rtag


def retr_gaus(gdat, gdatmodi, indxsamp, stdv, inpl=False):
    
    if gdat.probrand > 0.:
        if rand() < gdat.probrand:
            gdatmodi.nextsamp[indxsamp] = rand()
            thisbool = False
        else:
            thisbool = True
    else:
        thisbool = True
    
    if thisbool:
        if inpl:
            samp = gdatmodi.nextsamp
        else:
            samp = gdatmodi.thissamp
        if isinstance(stdv, float):
            gdatmodi.nextsamp[indxsamp] = samp[indxsamp] + normal(scale=stdv)
        else:
            for k in range(stdv.size):
                gdatmodi.nextsamp[indxsamp[k]] = samp[indxsamp[k]] + normal(scale=stdv[k])
            
       
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
    
    gang = arccos(cos(lgal) * cos(bgal))

    return gang


def retr_aang(lgal, bgal):

    aang = arctan2(bgal, lgal)

    return aang


def show_samp(gdat, gdatmodi):
    print '%22s %14s %14s %14s %14s %14s %14s %10s' % ('name', 'thissamp', 'nextsamp', 'thissampvarb', 'nextsampvarb', 'diffsampvarb', 'prop', 'indxstdp')
    for k in gdat.indxpara:
        if k == gdat.numbfixp:
            print
        if k < gdat.numbfixp:
            name = gdat.namefixp[k]
        else:
            name = ''

        gdatmodi.indxparaprop = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp)))
        if k in concatenate(gdatmodi.thisindxsamplgal):
            print

        try:
            strgboolmodi = '%s' % (k in gdatmodi.indxsampmodi)
        except:
            strgboolmodi = ''
        print '%22s %14.4g %14.4g %14.4g %14.4g %14.4g %14s %10d' % (name, getattr(gdatmodi, 'thissamp')[k], getattr(gdatmodi, 'nextsamp')[k], gdatmodi.thissampvarb[k], \
                                               gdatmodi.nextsampvarb[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k], strgboolmodi, gdatmodi.indxstdppara[k])
    
    for strg in ['this', 'next']:
        print '%sindxpntsfull' % strg
        print getattr(gdatmodi, strg + 'indxpntsfull')
        print '%sindxsamplgal' % strg
        print getattr(gdatmodi, strg + 'indxsamplgal')
        print '%sindxsampbgal' % strg
        print getattr(gdatmodi, strg + 'indxsampbgal')
        if gdat.numbener > 1:
            print '%sindxsampsind' % strg
            print getattr(gdatmodi, strg + 'indxsampsind')
            print '%sindxsampcurv' % strg
            print getattr(gdatmodi, strg + 'indxsampcurv')
            print '%sindxsampexpo' % strg
            print getattr(gdatmodi, strg + 'indxsampexpo')
        print '%sindxsampcomp' % strg
        print getattr(gdatmodi, strg + 'indxsampcomp')
        print


def retr_prop(gdat, gdatmodi):
 
    if gdat.verbtype > 1:
        print 'retr_prop()'
   
    for k in gdat.indxprop:
        retr_gaus(gdat, gdatmodi, gdat.indxfixpprop[k], gdatmodi.stdvstdp[k])
    
    for k in gdat.indxiact:
        gdatmodi.nextsampvarb[gdat.indxfixpiact[k]] = gdatmodi.thissampvarb[gdat.indxfixpiact[k]]
   
    for k in gdat.indxfixpdist:
        gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, '', gdatmodi.nextsamp[k], k)

    # rescale the unit sample vector if a hyperparameter controlling the distribution of PS properties is being updated
    for l in gdat.indxpopl: 
        
        ## flux distribution
        if gdat.fluxdisttype[l] == 'powr':
            fluxunit = cdfn_flux_powr(gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]], gdat.minmflux, gdat.maxmflux, gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslop[l]])
        if gdat.fluxdisttype[l] == 'brok':
            fluxunit = cdfn_flux_brok(gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]], gdat.minmflux, gdat.maxmflux, gdatmodi.nextsampvarb[gdat.indxfixpfluxdistbrek[l]], \
                                            gdatmodi.nextsampvarb[gdat.indxfixpfluxdistsloplowr[l]], gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslopuppr[l]])
        
        gdatmodi.nextsamp[gdatmodi.thisindxsampflux[l]] = fluxunit
    
        if gdat.numbener > 1:
            sindunit = cdfn_gaus(gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[l]], gdatmodi.nextsampvarb[gdat.indxfixpsinddistmean[l]], \
                                                                                                                gdatmodi.nextsampvarb[gdat.indxfixpsinddiststdv[l]])
            gdatmodi.nextsamp[gdatmodi.thisindxsampsind[l]] = sindunit
            
            if gdat.spectype[l] == 'curv':
                curvunit = cdfn_gaus(gdatmodi.thissampvarb[gdatmodi.thisindxsampcurv[l]], gdatmodi.nextsampvarb[gdat.indxfixpcurvdistmean[l]], \
                                                                                                                gdatmodi.nextsampvarb[gdat.indxfixpcurvdiststdv[l]])
                gdatmodi.nextsamp[gdatmodi.thisindxsampcurv[l]] = curvunit

            if gdat.spectype[l] == 'expo':
                expounit = cdfn_gaus(gdatmodi.thissampvarb[gdatmodi.thisindxsampexpo[l]], gdatmodi.nextsampvarb[gdat.indxfixpexpodistmean[l]], \
                                                                                                                gdatmodi.nextsampvarb[gdat.indxfixpexpodiststdv[l]])
                gdatmodi.nextsamp[gdatmodi.thisindxsampexpo[l]] = expounit


    # number of unit sample vector elements to be modified
    numbcompmodi = gdat.numbcomp[gdatmodi.indxpoplmodi]
    
    gdatmodi.nextindxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
    
    if gdatmodi.propbrth:
        
        # find an empty slot in the PS list
        for k in range(gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]):
            if not k in gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi]:
                indxpntsbrth = k
                break

        # sample indices to add the new PS
        gdatmodi.indxsamptran = gdat.indxsampcomp[0] + gdat.numbtrapcuml[gdatmodi.indxpoplmodi] + indxpntsbrth * gdat.numbcomp[gdatmodi.indxpoplmodi] + \
                                                                                                                                 gdat.indxcomp[gdatmodi.indxpoplmodi]
        
        # sample auxiliary variables
        gdatmodi.auxipara = rand(numbcompmodi)
        
        gdatmodi.nextsamp[gdatmodi.indxsamptran] = gdatmodi.auxipara

        if gdat.verbtype > 1:
            print 'auxipara'
            print gdatmodi.auxipara
            print 'numbcompmodi'
            print numbcompmodi
            print 'numbcompmodi'
            print numbcompmodi
                
    # death
    if gdatmodi.propdeth:
        
        # occupied PS index to be killed
        dethindxindxpnts = choice(arange(gdatmodi.thissamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        indxpntsdeth = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][dethindxindxpnts]

        # sample indices to add the new PS
        gdatmodi.indxsamptran = gdat.indxsampcomp[0] + gdat.numbtrapcuml[gdatmodi.indxpoplmodi] + indxpntsdeth * gdat.numbcomp[gdatmodi.indxpoplmodi] + \
                                                                                                                                 gdat.indxcomp[gdatmodi.indxpoplmodi]
        
        if gdat.verbtype > 1:
            print 'dethindxpnts: ', indxpntsdeth
            print 'dethindxindxpnts: ', dethindxindxpnts
            print
  
    ## birth
    if gdatmodi.propbrth or gdatmodi.propsplt:
        
        # change the number of PS
        gdatmodi.nextsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].append(indxpntsbrth)

    ## death
    if gdatmodi.propdeth or gdatmodi.propmerg:
        temp = deepcopy(gdatmodi.thisindxpntsfull)[gdatmodi.indxpoplmodi]
        indxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
        indxpntsfull[gdatmodi.indxpoplmodi].remove(indxpntsdeth)
        
        # change the number of PS
        gdatmodi.nextsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
        
        # remove the PS from the occupied PS list
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].remove(indxpntsdeth)
    else:
        indxpntsfull = gdatmodi.thisindxpntsfull
    
    gdatmodi.dictmodi['indxsamplgal'], gdatmodi.dictmodi['indxsampbgal'], gdatmodi.dictmodi['indxsampflux'], gdatmodi.dictmodi['indxsampsind'], \
        	    gdatmodi.dictmodi['indxsampcurv'], gdatmodi.dictmodi['indxsampexpo'], gdatmodi.dictmodi['indxsampcomp'] = retr_indx(gdat, indxpntsfull, gdat.spectype)
    
    # PSs
    gdatmodi.lfctprop = 0.
    for l in gdat.indxpopl:
        for k in range(len(indxpntsfull[l])):
            for f in gdat.indxcomp[l]:
                gdatmodi.thisstdv = gdatmodi.stdvstdp[gdat.indxstdpcomp[f]] / (gdatmodi.thissampvarb[gdatmodi.dictmodi['indxsampflux'][l][k]] / gdat.minmflux)**0.5
                gdatmodi.nextstdv = gdatmodi.stdvstdp[gdat.indxstdpcomp[f]] / (gdatmodi.nextsampvarb[gdatmodi.dictmodi['indxsampflux'][l][k]] / gdat.minmflux)**0.5
                retr_gaus(gdat, gdatmodi, gdatmodi.dictmodi['indxsamp' + gdat.liststrgcomp[l][f]][l][k], gdatmodi.thisstdv)
                gdatmodi.lfctprop += sum((gdatmodi.nextsamp[gdat.indxfixpprop] - gdatmodi.thissamp[gdat.indxfixpprop]) * \
                                                                                            (1. / gdatmodi.nextstdv**2 - 1. / gdatmodi.thisstdv**2))
    
    print 'gdatmodi.lfctprop'
    print gdatmodi.lfctprop
    gdatmodi.lfctprop = 0.

    if gdat.numbtrap > 0:
        gdatmodi.indxsampchec = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.dictmodi['indxsampcomp'])))
        if gdatmodi.propbrth:
            gdatmodi.indxsampmodi = concatenate((gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsampchec, gdatmodi.indxsamptran))
        else:
            gdatmodi.indxsampmodi = concatenate((gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsampchec))
    else:
        gdatmodi.indxsampchec = gdat.indxfixpprop
        gdatmodi.indxsampmodi = gdat.indxfixpprop
    
    # reject the sample if proposal is outside the prior
    indxchecfail = where((gdatmodi.nextsamp[gdatmodi.indxsampchec] < 0.) | (gdatmodi.nextsamp[gdatmodi.indxsampchec] > 1.))[0]
    if indxchecfail.size > 0:
        if gdat.verbtype > 1:
            print 'Proposal rejected due to proposal outside the prior'
            show_samp(gdat, gdatmodi)
            print 'indxchecfail'
            print indxchecfail
            print
        gdatmodi.boolreje = True
    
    if not gdatmodi.boolreje:
        #gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = retr_sampvarb(gdat, gdatmodi.nextindxpntsfull, gdatmodi.nextsamp, 'next')[gdatmodi.indxsampmodi]
        gdatmodi.nextsampvarb = retr_sampvarb(gdat, gdatmodi.nextindxpntsfull, gdatmodi.nextsamp, 'next')
   
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbpntsmodi = 3
        
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        gdatmodi.spltindxindxpnts = choice(thisindxindxpnts)
    
        # update the sample vector
        gdatmodi.indxsampfrst = gdat.indxsampcomp[0] + gdat.numbtrap * gdatmodi.indxpoplmodi + \
                                                int(gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]) * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp[gdatmodi.indxpoplmodi]

        gdatmodi.indxsampseco = gdat.indxsampcomp[0] + gdat.numbtrap * gdatmodi.indxpoplmodi + indxpntsbrth * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlseco = gdatmodi.indxsampseco + gdat.numbcomp[gdatmodi.indxpoplmodi]
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, indxfinlfrst, dtype=int), arange(gdatmodi.indxsampseco, indxfinlseco, dtype=int)))
        
        thislgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisbgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, gdatmodi.spltindxindxpnts]]
        gdatmodi.fluxpare = thisspec[gdat.indxenerfluxdist[0]]
        thisspep = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts, :]]
        
        # determine the new components
        # temp -- only valid for power-law energy spectrum
        gdatmodi.auxipara = empty(gdat.numbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.auxipara[0] = rand()
        gdatmodi.auxipara[1] = rand() * gdat.radispmr
        gdatmodi.auxipara[2] = rand() * 2. * pi
        # temp
        if gdat.numbener > 1:
            gdatmodi.auxipara[3] = icdf_gaus(rand(), gdat.sinddistmean[gdatmodi.indxpoplmodi], gdat.sinddiststdv[gdatmodi.indxpoplmodi])
            
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

        if fabs(gdatmodi.spltlgalfrst) > gdat.maxmgang or fabs(gdatmodi.spltlgalseco) > gdat.maxmgang or \
                                            fabs(gdatmodi.spltbgalfrst) > gdat.maxmgang or fabs(gdatmodi.spltbgalseco) > gdat.maxmgang or \
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
                raise Exception('Number of pairs should not be zero in the reverse proposal of a split')

            ## first new component
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcomplgal] = cdfn_self(gdatmodi.spltlgalfrst, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompbgal] = cdfn_self(gdatmodi.spltbgalfrst, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompflux] = cdfn_flux_powr(gdatmodi.fluxfrst, gdat.minmflux, gdat.maxmflux, \
                                                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompsind] = cdfn_gaus(gdatmodi.spltsindfrst, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
    
            nextspecfrst = retr_spec(gdat, gdatmodi.fluxfrst, spep=gdatmodi.spltsindfrst, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

            ## second new component
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcomplgal] = cdfn_self(gdatmodi.spltlgalseco, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcompbgal] = cdfn_self(gdatmodi.spltbgalseco, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcompflux] = cdfn_flux_powr(gdatmodi.fluxseco, gdat.minmflux, gdat.maxmflux, \
                                                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcompsind] = cdfn_gaus(gdatmodi.spltsindseco, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
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
            gdatmodi.auxipara = zeros(gdat.numbcomp[gdatmodi.indxpoplmodi])
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
            gdatmodi.nextsamp[gdatmodi.indxsampfrst] = cdfn_self(gdatmodi.lgalpare, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+1] = cdfn_self(gdatmodi.bgalpare, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+2] = cdfn_flux_powr(gdatmodi.fluxpare, gdat.minmflux, gdat.maxmflux, \
                                                                                        gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+3] = gdatmodi.thissamp[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][mergindxindxpntsfrst, :]]

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

    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and not gdatmodi.boolreje:
        
        ## Jacobian
        jcbnfacttemp = log(gdatmodi.fluxpare * fabs(gdatmodi.auxipara[1] * (sin(gdatmodi.auxipara[2]) * cos(gdatmodi.auxipara[2]) + cos(gdatmodi.auxipara[2])**2)))
        if gdatmodi.propsplt:
            gdatmodi.jcbnfact = jcbnfacttemp
        else:
            gdatmodi.jcbnfact = -jcbnfacttemp
        
        ## combinatorial factor
        thisnumbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]
        if gdatmodi.propsplt:
            gdatmodi.combfact = log(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2 / gdatmodi.nextnumbpair)
        else:
            gdatmodi.combfact = log(gdatmodi.thisnumbpair / gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2)
        
        if gdat.verbtype > 1:
            print 'jcbnfact'
            print gdatmodi.jcbnfact
            print 'combfact'
            print gdatmodi.combfact
            print

    else:
        gdatmodi.jcbnfact = 0.
        gdatmodi.combfact = 0.
   
        
def retr_factoaxi(gdat, bins, norm, indx):

    factoaxi = 1. + norm[:, None, None] * (bins[None, None, :] / gdat.oaxipivt)**indx[:, None, None]
    
    return factoaxi


def gang_detr():

    gang, aang, lgal, bgal = sympy.symbols('gang aang lgal bgal')

    AB = sympy.matrices.Matrix([[a1*b1,a1*b2,a1*b3],[a2*b1,a2*b2,a2*b3],[a3*b1,a3*b2,a3*b3]])
    print AB.det()


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, binsoaxi=None, varioaxi=None, strgpara=''):

    numbpsfpform = getattr(gdat, strgpara + 'numbpsfpform')
    numbpsfptotl = getattr(gdat, strgpara + 'numbpsfptotl')
    
    indxpsfpinit = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    if varioaxi:
        indxpsfponor = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp]
        indxpsfpoind = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp] + 1

    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * arcsin(sqrt(2. - 2. * cos(gdat.binsangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        if varioaxi:
            scalangl = thisangl[None, :, None, None]
        else:
            scalangl = thisangl[None, :, None]
    
    if varioaxi:
        factoaxi = retr_factoaxi(gdat, binsoaxi, psfp[indxpsfponor], psfp[indxpsfpoind])
   
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
        indxwvec = where((gdat.meanwvec > gdat.binswvecodimplot[k]) & (gdat.meanwvec < gdat.binswvecodimplot[k+1]))
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


def retr_varb(gdat, gdatmodi, strg, perc='medi'):
        
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
                    
                                ## compute distance

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

        if attr == 'thissampvarb' or attr == 'thissamp':
            
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
        cnts[:, k, :] += spec[:, k, None] * gdat.expo[:, indxpixltemp[k], :]
    if gdat.enerdiff:
        cnts *= gdat.deltener[:, None, None]
    
    return cnts


def setpinit(gdat, boolinitsetp=False):

    # samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    
    # samples to be saved from all chains
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
        gdat.pathpostlpri = gdat.pathpost + 'lpri/'
        gdat.pathpostfixp = gdat.pathpost + 'fixp/'
        gdat.pathpostfixpproc = gdat.pathpostfixp + 'proc/'
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
    ## make the directories 
    for attr, valu in gdat.__dict__.iteritems():
        if attr.startswith('path'):
            os.system('mkdir -p %s' % valu)
 
    # plotting
    ## number of bins in histogram plots
    gdat.numbbinsplot = 20
    ## number of bins in hyperprior plots
    gdat.numbbinsplotprio = 100
    # temp
    ### number of standard deviations away from the mean of spectral parameters to plot
    gdat.numbstdvspepdist = 3.
    ### minima and maxima for spectral parameters
    gdat.numbstdv = 3.

    # list of source components, i.e., primary model parameters
    gdat.liststrgcompdefa = ['lgal', 'bgal', 'flux', 'sind', 'curv', 'expo']
    
    # total maximum number of PS
    gdat.maxmnumbpntstotl = sum(gdat.maxmnumbpnts)
    gdat.indxpntstotl = arange(gdat.maxmnumbpntstotl)
    gdat.maxmnumbpntscumr = cumsum(gdat.maxmnumbpnts)
    gdat.maxmnumbpntscuml = concatenate((array([0]), gdat.maxmnumbpntscumr[:-1]))
   
    gdat.lablener = 'E'
    gdat.lablenertotl = '$%s$ [%s]' % (gdat.lablener, gdat.strgenerunit)
    gdat.lablgang = r'\theta'
    gdat.lablaang = r'\phi'
    gdat.labllgalunit = gdat.lablgangunit
    gdat.lablbgalunit = gdat.lablgangunit
    gdat.lablaangunit = ''
    gdat.lablsind = 's'
    gdat.lablsindunit = ''
    gdat.lablcurv = r'\kappa'
    gdat.lablcurvunit = ''
    gdat.lablexpo = r'\epsilon'
    gdat.lablexpounit = gdat.strgenerunit
    gdat.lablspec = gdat.lablflux
    gdat.lablspecunit = gdat.lablfluxunit
    gdat.lablcnts = 'C'
    gdat.lablcntsunit = ''
    
    # set model sample vector indices
    retr_indxsamp(gdat)
    
    gdat.lablfeat = {}
    gdat.dictglob = {}
    gdat.lablfeatunit = {}
    gdat.lablfeattotl = {}
    for strgfeat in gdat.liststrgfeat:
        gdat.lablfeat[strgfeat] = getattr(gdat, 'labl' + strgfeat)
        gdat.lablfeatunit[strgfeat] = getattr(gdat, 'labl' + strgfeat + 'unit')
        if gdat.lablfeatunit[strgfeat] != '':
            gdat.lablfeatunit[strgfeat] = ' [%s]' % gdat.lablfeatunit[strgfeat]
        gdat.lablfeattotl[strgfeat] = '$%s$%s' % (gdat.lablfeat[strgfeat], gdat.lablfeatunit[strgfeat])
        if strgfeat == 'flux' and gdat.pntstype == 'lens' or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal':
            gdat.dictglob['fact' + strgfeat + 'plot'] = gdat.anglfact
        else:
            gdat.dictglob['fact' + strgfeat + 'plot'] = 1.
        setattr(gdat, 'numb' + strgfeat + 'plot', 20)
        
        if strgfeat == 'flux':
            gdat.dictglob['scal' + strgfeat + 'plot'] = 'logt'
        else:
            gdat.dictglob['scal' + strgfeat + 'plot'] = 'self'

    # sweeps to be saved
    gdat.indxswep = arange(gdat.numbswep)
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
        gdat.numblpri = 1 + 7 * gdat.numbpopl
    else:
        gdat.numblpri = 0

    # size of the auxiliary variable propobability density vector
    gdat.numblpau = gdat.maxmnumbcomp * gdat.numbpopl
    
    # process index
    gdat.indxproc = arange(gdat.numbproc)

    # flag to indicate whether information from a deterministic catalog will be used or not
    # temp -- if datatype == 'inpt' trueinfo should depend on whether truexxxx are provided
    gdat.trueinfo = gdat.datatype == 'mock' or gdat.exprinfo
    
    if gdat.pntstype == 'lens':
        gdat.minmmass = retr_massfrombein(gdat.minmflux)
        gdat.maxmmass = retr_massfrombein(gdat.maxmflux)
        retr_axis(gdat, 'mass', gdat.minmmass, gdat.maxmmass, gdat.numbbinsplot)
        retr_axis(gdat, 'bein', gdat.minmflux, gdat.maxmflux, gdat.numbbinsplot)

    gdat.numbsinddistpara = 2
    gdat.numbfluxdistpara = 4
    
    gdat.listlablcompfrac = deepcopy(gdat.nameback)
    if gdat.pntstype == 'lght':
        gdat.listlablcompfrac.append('PS')
    if gdat.pntstype == 'lens':
        gdat.listlablcompfrac.append('Source')
        gdat.listlablcompfrac.append('Host')
    
    gdat.listlablcompfracspec = deepcopy(gdat.listlablcompfrac)
    gdat.listlablcompfracspec += ['Data']
    if len(gdat.listlablcompfrac) > 1:
        gdat.listlablcompfracspec.append('Total Model')
    
    gdat.numblablcompfrac = len(gdat.listlablcompfrac)
    gdat.numblablcompfracspec = len(gdat.listlablcompfracspec)
    
    if gdat.numbener > 1:
        gdat.enerfluxdist = gdat.meanener[gdat.indxenerfluxdist]
        if gdat.enerfluxdist == 0.:
            raise Exception('Pivot energy cannot be zero.')
        gdat.factspecener = (gdat.meanener / gdat.enerfluxdist)**(-sqrt(amin(gdat.minmsinddistmean) * amax(gdat.maxmsinddistmean)))
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
    gdat.offstext = gdat.maxmgangdata * 0.05
    ## figure size
    gdat.plotsize = 7
    ## size of the images
    gdat.sizeimag = 1.3 * gdat.plotsize
    ## text
    if gdat.datatype == 'mock':
        gdat.truelabl = 'Mock'
    if gdat.datatype == 'inpt':
        gdat.truelabl = gdat.strgcatl

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
        gdat.maxmoaxi = 1.1 * sqrt(2.) * gdat.maxmgang
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

    # the function to measure time
    # temp
    gdat.strgfunctime = 'clck'
    if gdat.strgfunctime == 'clck':
        gdat.functime = time.clock
    if gdat.strgfunctime == 'time':
        gdat.functime = time.time

    # axes
    gdat.minmlgaldata = -gdat.maxmgangdata
    gdat.maxmlgaldata = gdat.maxmgangdata
    gdat.minmbgaldata = -gdat.maxmgangdata
    gdat.maxmbgaldata = gdat.maxmgangdata
    ## longitude
    gdat.numblgalpntsprob = gdat.numbsidepntsprob
    gdat.numbbgalpntsprob = gdat.numbsidepntsprob
    gdat.binslgalpntsprob = linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.binsbgalpntsprob = linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.indxlgalpntsprob = arange(gdat.numblgalpntsprob)
    gdat.indxbgalpntsprob = arange(gdat.numbbgalpntsprob)

    if gdat.pntstype == 'lens': 
        retr_axis(gdat, 'deflplot', -gdat.maxmgang, gdat.maxmgang, 50)

    # lensing problem setup
    ## number of subhalos to plot
    gdat.numbdeflpnts = 3
    ## number of deflection components to plot
    gdat.numbdeflsing = gdat.numbdeflpnts + 2

    # convenience variables
    gdat.indxfluxprox = arange(gdat.numbfluxprox)
    gdat.binsfluxprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
    gdat.meanfluxprox = sqrt(gdat.binsfluxprox[1:] * gdat.binsfluxprox[:-1])
    
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
        gdat.apix = (2. * gdat.maxmgangdata)**2
    else:
        if gdat.pixltype == 'cart':
            gdat.binslgalcart = linspace(gdat.minmlgaldata, gdat.maxmlgaldata, gdat.numbsidecart + 1)
            gdat.binsbgalcart = linspace(gdat.minmbgaldata, gdat.maxmbgaldata, gdat.numbsidecart + 1)
            gdat.lgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
            gdat.bgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
            gdat.apix = (2. * gdat.maxmgangdata / gdat.numbsidecart)**2
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
   
            gdat.indxpixlrofi = where((fabs(lgalheal) < gdat.maxmgangdata) & (fabs(bgalheal) < gdat.maxmgangdata))[0]
            
            gdat.indxpixlrofimarg = where((fabs(lgalheal) < 1.2 * gdat.maxmgang) & (fabs(bgalheal) < 1.2 * gdat.maxmgang))[0]

            gdat.lgalgrid = lgalheal
            gdat.bgalgrid = bgalheal

        gdat.indxpixlfull = arange(gdat.numbpixlfull)

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
        
        retr_axis(gdat, 'wvecodimplot', gdat.minmwvecodim, gdat.maxmwvecodim, gdat.numbwvecodim)
         
        gdat.numbsidewvec = gdat.numbsidecart / 2
        temp = fft.fftfreq(gdat.numbsidewvec, gdat.sizepixl)
        gdat.meanwveclgal, gdat.meanwvecbgal = meshgrid(temp, temp, indexing='ij')
        gdat.meanwveclgal *= gdat.factkpcs
        gdat.meanwvecbgal *= gdat.factkpcs
        gdat.meanwvec = sqrt(gdat.meanwveclgal**2 + gdat.meanwvecbgal**2)

    # element parameter vector indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompflux = 2
    gdat.indxcompsind = 3
    gdat.indxcompcurv = 4
    gdat.indxcompexpo = 4
    gdat.indxcomp = [[] for l in gdat.indxpopl]
    for l in gdat.indxpopl:
        gdat.indxcomp[l] = arange(gdat.numbcomp[l])
    gdat.indxpnts = []
    for l in gdat.indxpopl:
        gdat.indxpnts.append(arange(gdat.maxmnumbpnts[l]))

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
    # temp
    if False and gdat.datatype == 'inpt':
        gdat.exprdataflux = gdat.exprdataflux[gdat.indxcubeincl]

        # find blobs in the given image
        mapstemp = tdpy.util.retr_cart(gdat.exprdataflux[0, :, 0], indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                            minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                            minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata)
        gdat.listblob = blob_doh(mapstemp, max_sigma=30, threshold=.01)
        for blob in gdat.listblob:
            y, x, r = blob
            print 'x'
            print x
            print 'y'
            print y
            print 'r'
            print r
            print

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
        path = gdat.pathpixlcnvt + 'pixlcnvt_%09g.p' % gdat.maxmgangdata

        if os.path.isfile(path):
            if gdat.verbtype > 0 and boolinitsetp:
                print 'Reading %s...' % path
            fobj = open(path, 'rb')
            gdat.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            gdat.pixlcnvt = zeros(gdat.numbpixlfull, dtype=int) - 1

            numbpixlmarg = gdat.indxpixlrofimarg.size
            for k in range(numbpixlmarg):
                dist = retr_angldistunit(gdat, lgalheal[gdat.indxpixlrofimarg[k]], bgalheal[gdat.indxpixlrofimarg[k]], gdat.indxpixl)
                gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)
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
                gdat.maxmangleval[h] = 3. * gdat.maxmgang
            else:   
                frac = gdat.specfraceval * gdat.binsfluxprox[0] / gdat.binsfluxprox[h+1]
                psfnwdth = retr_psfnwdth(gdat, gdat.exprpsfn, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathpixlprox + 'indxpixlprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), \
                                                                                                            1e2 * amax(gdat.maxmangleval), gdat.numbfluxprox)
        if gdat.verbtype > 0 and boolinitsetp:
            print 'PSF evaluation will be performed up to %.3g %s for the largest flux.' % (amax(gdat.maxmangleval) * gdat.anglfact, gdat.lablfeatunit['gang'])
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
    gdat.limtpntshist = [0.5, 10**ceil(log10(gdat.maxmnumbpntstotl))]
    
    ## marker size
    gdat.minmmrkrsize = 50
    gdat.maxmmrkrsize = 500
    
    ## ROI
    gdat.exttrofi = array([gdat.minmlgaldata, gdat.maxmlgaldata, gdat.minmbgaldata, gdat.maxmbgaldata])
    gdat.exttrofi *= gdat.anglfact 
    gdat.frambndrdata = gdat.maxmgangdata * gdat.anglfact
    gdat.frambndrmodl = gdat.maxmgang * gdat.anglfact

    ## marker opacity
    gdat.alphmrkr = 0.5
    gdat.alphpnts = 0.4
    gdat.alphmaps = 1.
    
    # number of colorbar ticks in the maps
    gdat.numbtickcbar = 11
    
    gdat.minmconv = 1e-2
    gdat.maxmconv = 1e1
    gdat.minmdeflcomp = 0.
    gdat.maxmdeflcomp = 0.1
    
    if gdat.numbener > 1:
        # temp
        gdat.minmexpo = gdat.minmexpodistmean
        gdat.maxmexpo = gdat.maxmexpodistmean
    
    ## plot limits for element parameters
    for strgcomp in gdat.liststrgcomptotl:
        for strglimt in ['minm', 'maxm']:
            try:
                truelimt = getattr(gdat, 'true' + strglimt + strgcomp)
            except:
                truelimt = None
             
            if strgcomp in ['sind', 'curv']:
                if strglimt == 'minm':
                    limt = getattr(gdat, 'minm' + strgcomp + 'distmean') - getattr(gdat, 'maxm' + strgcomp + 'diststdv')
                else:
                    limt = getattr(gdat, 'maxm' + strgcomp + 'distmean') + getattr(gdat, 'maxm' + strgcomp + 'diststdv')
            else:
                limt = getattr(gdat, strglimt + strgcomp)
            
            setattr(gdat, strglimt + strgcomp, limt)

            if truelimt == None:
                setattr(gdat, strglimt + strgcomp + 'plot', limt)
            else:
                limt = min(limt, truelimt)
                setattr(gdat, strglimt + strgcomp + 'plot', limt)
    
    gdat.minmaang = -pi
    gdat.maxmaang = pi

    # temp
    gdat.minmaangplot = -pi
    gdat.maxmaangplot = pi
    
    gdat.minmgangplot = 0.
    gdat.maxmgangplot = gdat.maxmlgal
    gdat.minmspecplot = gdat.minmfluxplot * gdat.factspecener
    gdat.maxmspecplot = gdat.maxmfluxplot * gdat.factspecener
    gdat.minmcntsplot = 0.1 * gdat.minmflux * mean(mean(gdat.expo, 1), 1)
    gdat.maxmcntsplot = gdat.maxmflux * mean(mean(gdat.expo, 1), 1)
    if gdat.enerdiff:
        gdat.minmcntsplot *= gdat.deltener * gdat.factspecener
        gdat.maxmcntsplot *= gdat.deltener * gdat.factspecener

    if gdat.pntstype == 'lens':
        liststrgcbar = ['conv', 'deflcomp']
        for strgcbar in liststrgcbar:
            retr_ticklabl(gdat, strgcbar)
   
    # temp
    gdat.minmspec = gdat.minmflux * gdat.factspecener
    gdat.maxmspec = gdat.maxmflux * gdat.factspecener

    for strgfeat in gdat.liststrgfeat:
        if strgfeat == 'spec':
            gdat.minmspecplot = gdat.minmfluxplot * gdat.factspecener
            gdat.maxmspecplot = gdat.maxmfluxplot * gdat.factspecener
            gdat.binsspecplot = gdat.binsfluxplot[None, :] * gdat.factspecener[:, None]
            gdat.meanspecplot = empty((gdat.numbener, gdat.numbbinsplot))
            for i in gdat.indxener:
                gdat.meanspecplot[i, :] = sqrt(gdat.binsspecplot[i, 1:] * gdat.binsspecplot[i, :-1])
        elif strgfeat == 'cnts':
            gdat.minmcnts = 0.1 * gdat.minmflux * mean(mean(gdat.expo, 1), 1)
            gdat.maxmcnts = gdat.maxmflux * mean(mean(gdat.expo, 1), 1)
            if gdat.enerdiff:
                gdat.minmcnts *= gdat.deltener * gdat.factspecener
                gdat.maxmcnts *= gdat.deltener * gdat.factspecener
            retr_axis(gdat, 'cntspivt', gdat.minmcnts[gdat.indxenerfluxdist[0]], gdat.maxmcnts[gdat.indxenerfluxdist[0]], gdat.numbbinsplot, scal='logt')
            gdat.binscnts = empty((gdat.numbener, gdat.numbbinsplot + 1))
            for i in gdat.indxener:
                gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbbinsplot + 1) # [1]
            gdat.meancnts = sqrt(gdat.binscnts[:, :-1] * gdat.binscnts[:, 1:]) 
        else:
            if strgfeat in ['flux', 'expo']:
                scal = 'logt'
            else:
                scal = 'self'
            maxm = getattr(gdat, 'maxm' + strgfeat + 'plot')
            minm = getattr(gdat, 'minm' + strgfeat + 'plot')
            retr_axis(gdat, strgfeat + 'plot', minm, maxm, gdat.numbbinsplot, scal=scal)
            retr_axis(gdat, strgfeat + 'plotprio', minm, maxm, gdat.numbbinsplotprio, scal=scal)
    
        gdat.dictglob['limt' + strgfeat + 'plot'] = [getattr(gdat, 'minm' + strgfeat), getattr(gdat, 'maxm' + strgfeat)]


def setpfinl(gdat, boolinitsetp=False):

    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = retr_cntsmaps(gdat, gdat.exprdataflux)
    
    if gdat.pixltype != 'unbd':
        gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix
        if gdat.enerdiff:
            gdat.datafluxmean /= gdat.deltener
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
       
    # spatially averaged background flux 
    gdat.backfluxmean = zeros((gdat.numbback, gdat.numbener))
    for c in gdat.indxback:
        for i in gdat.indxener:
            if gdat.correxpo:
                gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, :] * gdat.expo[i, :, :]) / sum(gdat.expo[i, :, :])
            else:
                gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, 0])
    
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
                    #if gdat.truenumbpopl != gdat.numbpopl or gdat.namefixp[k].startswith('fluxdist') and gdat.truefluxdisttype[l] != gdat.fluxdisttype[l]:
                    #    continue
                    strg = gdat.namefixp[k][:-4]
                    print 'true' + strg
                    gdat.corrfixp[k] = getattr(gdat, 'mock' + strg)[l]
                else:
                    strg = gdat.namefixp[k]
                    print 'true' + strg
                    gdat.corrfixp[k] = getattr(gdat, 'mock' + strg)
    
    if gdat.pixltype == 'unbd':
        gdat.bgalgrid = gdat.datacnts[0, :, 0, 0]
        gdat.lgalgrid = gdat.datacnts[0, :, 0, 1]

    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(10000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)
    gdat.indxcubesave = meshgrid(gdat.indxener, gdat.indxpixlsave, gdat.indxevtt, indexing='ij')
   
    # temp -- needed for experimental PSF
    # if gdat.evalpsfnpnts:
    #     gdat.truenumbpsfpform, gdat.truenumbpsfpoaxi, gdat.truenumbpsfptotl, gdat.trueindxpsfponor, gdat.trueindxpsfpoind = \
    #                                                                                         retr_indxpsfp(gdat, gdat.truepsfntype, gdat.truevarioaxi)
    #     if gdat.truevarioaxi:
    #         gdat.truefactoaxi = retr_factoaxi(gdat, gdat.binsoaxi, gdat.truepsfp[gdat.trueindxpsfponor], gdat.truepsfp[gdat.trueindxpsfpoind])
    
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
            for strgtype in gdat.liststrgtype:
                print getattr(gdat, strgpara + strgtype)
            if strgpara == '':
                listfeat = ('name', 'strg', 'scal', 'minm', 'maxm', 'true')
                print '%20s%25s%5s%20s%20s%20s' % listfeat
            else:
                listfeat = ('name', 'strg', 'scal', 'minm', 'maxm', 'sampvarb')
                print '%20s%25s%5s%20s%20s%20s' % listfeat
            for k in getattr(gdat, strgpara + 'indxfixp'):
                if strgpara == '':
                    print '%20s%25s%5s%20.6g%20.6g' % (gdat.namefixp[k], gdat.strgfixp[k], gdat.scalfixp[k], gdat.minmfixp[k], gdat.maxmfixp[k])
                else:
                    print '%20s%25s%5s%20.6g%20.6g%20.6g' % (gdat.namefixp[k], gdat.strgfixp[k], gdat.scalfixp[k], gdat.minmfixp[k], gdat.maxmfixp[k], gdat.truefixp[k])
    
    if gdat.trueinfo and gdat.correxpo and gdat.pntstype == 'lght':
        gdat.exprfwhm = 2. * retr_psfnwdth(gdat, gdat.exprpsfn, 0.5)

    # proposals
    gdat.indxfixpprop = []
    for k, strg in enumerate(gdat.namefixp):
        if k in gdat.indxfixpnumbpnts:
            thisbool = False
        else:
            if k in gdat.indxfixphypr:
                strgtemp = 'hypr'
            elif k in gdat.indxfixppsfp:
                strgtemp = 'psfp'
            elif k in gdat.indxfixpbacp:
                strgtemp = 'bacp'
            elif k in gdat.indxfixplenp:
                strgtemp = 'lenp'
            thisbool = getattr(gdat, 'prop' + strgtemp)
            
        if thisbool:
            gdat.indxfixpprop.append(getattr(gdat, 'indxfixp' + strg))
    gdat.indxfixpprop = array(gdat.indxfixpprop) 
    gdat.numbfixpprop = gdat.indxfixpprop.size
    gdat.indxprop = arange(gdat.numbfixpprop)
    
    # parameters not subject to proposals
    gdat.indxfixpiact = setdiff1d(gdat.indxfixp, gdat.indxfixpprop)
    gdat.numbfixpiact = gdat.indxfixpiact.size
    gdat.indxiact = arange(gdat.numbfixpiact)
    
    gdat.strgproptype = array([])
    gdat.nameproptype = array([])
   
    if gdat.probtran == None:
        if gdat.numbtrap > 0:
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
       
    cntr = tdpy.util.cntr()
    if gdat.numbtrap > 0.:
    
        if gdat.probtran > 0.:
            # birth
            gdat.indxproptypebrth = cntr.incr()
            gdat.strgproptype = append(gdat.strgproptype, r'$\mathcal{B}$')
            gdat.nameproptype = append(gdat.nameproptype, 'brth')
            
            # death
            gdat.indxproptypedeth = cntr.incr()
            gdat.strgproptype = append(gdat.strgproptype, r'$\mathcal{D}$')
            gdat.nameproptype = append(gdat.nameproptype, 'deth')
            
            if gdat.probbrde < 1.:
                # split
                gdat.indxproptypesplt = cntr.incr()
                gdat.strgproptype = append(gdat.strgproptype, r'$\mathcal{S}$')
                gdat.nameproptype = append(gdat.nameproptype, 'splt')
                
                # merge
                gdat.indxproptypemerg = cntr.incr()
                gdat.strgproptype = append(gdat.strgproptype, r'$\mathcal{M}$')
                gdat.nameproptype = append(gdat.nameproptype, 'merg')
    
    gdat.indxstdplgal = gdat.numbfixpprop
    gdat.indxstdpbgal = gdat.numbfixpprop + 1
    gdat.indxstdpflux = gdat.numbfixpprop + 2
    if gdat.numbener > 1:
        gdat.indxstdpsind = gdat.numbfixpprop + 3
        gdat.indxstdpcurv = gdat.numbfixpprop + 4
        gdat.indxstdpexpo = gdat.numbfixpprop + 4
    gdat.numbstdp = gdat.numbfixpprop + gdat.maxmnumbcomp
    gdat.strgstdp = concatenate((gdat.strgfixp[gdat.indxfixpprop], gdat.liststrgcomptotl))
    gdat.namestdp = concatenate((gdat.namefixp[gdat.indxfixpprop], gdat.liststrgcomptotl))
    gdat.indxstdp = arange(gdat.numbstdp)
    gdat.indxstdpfixp = arange(gdat.numbfixpprop)
    gdat.indxstdpcomp = setdiff1d(gdat.indxstdp, gdat.indxstdpfixp)
    
    # proposal scale
    gdat.stdvstdp = 1e-5 + zeros(gdat.numbstdp)
    
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
            #labl[k] = tdpy.util.mexp(sinh(tick[k]))
            labl[k] = '%.3g' % sinh(tick[k])
        else:
            #labl[k] = tdpy.util.mexp(tick[k])
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


def retr_indxsamp(gdat, strgpara=''):
    
    spectype = getattr(gdat, strgpara + 'spectype')
    spatdisttype = getattr(gdat, strgpara + 'spatdisttype')
    fluxdisttype = getattr(gdat, strgpara + 'fluxdisttype')
    varioaxi = getattr(gdat, strgpara + 'varioaxi')
    psfntype = getattr(gdat, strgpara + 'psfntype')

    numbback = getattr(gdat, strgpara + 'numbback')
    indxback = arange(numbback)
    
    numbbacp = numbback * gdat.numbener

    # population index vector
    numbpopl = getattr(gdat, strgpara + 'numbpopl')
    indxpopl = arange(numbpopl, dtype=int) 

    liststrgcomp = [[] for l in indxpopl]
    liststrgfeat = ['gang', 'aang', 'lgal', 'bgal', 'flux', 'spec', 'cnts']
    liststrgfeatprio = [[] for l in indxpopl]
    if gdat.numbener > 1 and gdat.pntstype == 'lght':
        liststrgfeat += ['sind']
    for l in indxpopl:
        liststrgcomp[l] = ['lgal', 'bgal', 'flux']
        if gdat.numbener > 1 and gdat.pntstype == 'lght':
            liststrgcomp[l] += ['sind']
        if spatdisttype[l] == 'gang':
            liststrgfeatprio[l] += ['gang', 'aang']
        elif spatdisttype[l] == 'disc' or spatdisttype[l] == 'unif':
            liststrgfeatprio[l] += ['lgal', 'bgal']
        liststrgfeatprio[l] += ['flux']
        if gdat.numbener > 1 and gdat.pntstype == 'lght':
            liststrgfeatprio[l] += ['sind']
            if spectype[l] == 'curv':
                liststrgcomp[l] += ['curv']
                liststrgfeatprio[l] += ['curv']
                if not 'curv' in liststrgfeat:
                    liststrgfeat += ['curv']
            if spectype[l] == 'expo':
                liststrgcomp[l] += ['expo']
                liststrgfeatprio[l] += ['expo']
                if not 'expo' in liststrgfeat:
                    liststrgfeat += ['expo']
    
    liststrgcomptotl = []
    for listsubb in liststrgcomp:
        for strg in listsubb:
            if not strg in liststrgcomptotl:
                liststrgcomptotl.append(strg)

    cntr = tdpy.util.cntr()
    
    dicttemp = {}
    if gdat.maxmnumbpntstotl > 0:
        for l in indxpopl:
            dicttemp['indxfixpnumbpntspop%d' % l] = cntr.incr()
        for l in indxpopl:
            dicttemp['indxfixpmeanpntspop%d' % l] = cntr.incr()
    
        liststrgvarb = ['gangdistscal', 'bgaldistscal', 'fluxdistslop', 'fluxdistbrek', 'fluxdistsloplowr', 'fluxdistslopuppr', 'sinddistmean', 'sinddiststdv', \
                                                                                                    'curvdistmean', 'curvdiststdv', 'expodistmean', 'expodiststdv']
        for strgvarb in liststrgvarb:
            strgtemp = 'indxfixp' + strgvarb
            dicttemp[strgtemp] = zeros(numbpopl, dtype=int) - 1
        
        for l in range(numbpopl):
            for strg in liststrgfeatprio[l]:
                if strg == 'gang' and spatdisttype[l] == 'gang':
                    dicttemp['indxfixpgangdistscalpop%d' % l] = cntr.incr()
                    dicttemp['indxfixpgangdistscal'][l] = dicttemp['indxfixpgangdistscalpop%d' % l]
                if strg == 'bgal' and spatdisttype[l] == 'disc':
                    dicttemp['indxfixpbgaldistscalpop%d' % l] = cntr.incr()
                    dicttemp['indxfixpbgaldistscal'][l] = dicttemp['indxfixpbgaldistscalpop%d' % l]
                if strg == 'flux':
                    if fluxdisttype[l] == 'powr':
                        dicttemp['indxfixpfluxdistsloppop%d' % l] = cntr.incr()
                        dicttemp['indxfixpfluxdistslop'][l] = dicttemp['indxfixpfluxdistsloppop%d' % l]
                    if fluxdisttype[l] == 'brok':
                        dicttemp['indxfixpfluxdistbrekpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpfluxdistsloplowrpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpfluxdistslopupprpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpfluxdistbrek'][l] = dicttemp['indxfixpfluxdistbrekpop%d' % l]
                        dicttemp['indxfixpfluxdistsloplowr'][l] = dicttemp['indxfixpfluxdistsloplowrpop%d' % l]
                        dicttemp['indxfixpfluxdistslopuppr'][l] = dicttemp['indxfixpfluxdistslopupprpop%d' % l]
                if gdat.numbener > 1:
                    if strg == 'sind':
                        dicttemp['indxfixpsinddistmeanpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpsinddiststdvpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpsinddistmean'][l] = dicttemp['indxfixpsinddistmeanpop%d' % l]
                        dicttemp['indxfixpsinddiststdv'][l] = dicttemp['indxfixpsinddiststdvpop%d' % l]
                    if strg == 'curv' and spectype[l] == 'curv':
                        dicttemp['indxfixpcurvdistmeanpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpcurvdiststdvpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpcurvdistmean'][l] = dicttemp['indxfixpcurvdistmeanpop%d' % l]
                        dicttemp['indxfixpcurvdiststdv'][l] = dicttemp['indxfixpcurvdiststdvpop%d' % l]
                    if strg == 'expo' and spectype[l] == 'expo':
                        dicttemp['indxfixpexpodistmeanpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpexpodiststdvpop%d' % l] = cntr.incr()
                        dicttemp['indxfixpexpodistmean'][l] = dicttemp['indxfixpexpodistmeanpop%d' % l]
                        dicttemp['indxfixpexpodiststdv'][l] = dicttemp['indxfixpexpodiststdvpop%d' % l]
    
    dicttemp['indxfixpnumbpnts'] = []
    dicttemp['indxfixpmeanpnts'] = []
    for strg, valu in dicttemp.iteritems():
        if strg[8:].startswith('numbpntsp'):
            dicttemp['indxfixpnumbpnts'].append(valu)
        if strg[8:].startswith('meanpntsp'):
            dicttemp['indxfixpmeanpnts'].append(valu)
    dicttemp['indxfixpnumbpnts'] = array(dicttemp['indxfixpnumbpnts'])

    dicttemp['indxfixpdist'] = []
    for strg, valu in dicttemp.iteritems():
        if strg[12:16] == 'dist' and isscalar(valu):
            dicttemp['indxfixpdist'].append(valu)
        
    dicttemp['indxfixphypr'] = dicttemp['indxfixpdist'] +  dicttemp['indxfixpmeanpnts']
    
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if psfntype == 'singgaus':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'singking':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'doubgaus':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'gausking':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'doubking':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
            if varioaxi:
                dicttemp['indxfixponorene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpoindene%devt%d' % (i, m)] = cntr.incr()
    
    dicttemp['indxfixponor'] = []
    dicttemp['indxfixpoind'] = []
    dicttemp['indxfixppsfp'] = []
    for strg, valu in dicttemp.iteritems():
        if strg.startswith('indxfixponore'):
            dicttemp['indxfixponor'].append(valu)
        if strg.startswith('indxfixpoinde'):
            dicttemp['indxfixpoind'].append(valu)
        if strg.startswith('indxfixpsigce') or strg.startswith('indxfixpsigte') or strg.startswith('indxfixpgamce') or strg.startswith('indxfixpgamte') or \
                                                                    strg.startswith('indxfixppsffe') or strg.startswith('indxfixponore') or strg.startswith('indxfixpoinde'):
            dicttemp['indxfixppsfp'].append(valu)
    dicttemp['indxfixppsfp'] = array(dicttemp['indxfixppsfp']) 
    numbpsfpform, numbpsfpoaxi, numbpsfptotl, indxpsfponor, indxpsfpoind = retr_indxpsfp(gdat, psfntype, varioaxi)
    
    numbpsfptotlevtt = gdat.numbevtt * numbpsfptotl
    numbpsfptotlener = gdat.numbener * numbpsfptotl
    numbpsfp = numbpsfptotl * gdat.numbener * gdat.numbevtt
    indxpsfpoaxi = arange(numbpsfpoaxi) 
    indxpsfpform = arange(numbpsfpform)
    indxpsfptotl = arange(numbpsfptotl)
   
    if varioaxi:
        indxfixppsfponor = indxfixppsfp[0] + indxpsfponor
        indxfixppsfpoind = indxfixppsfp[0] + indxpsfpoind
        indxfixppsfpoaxi = sort(concatenate((indxfixppsfponor, indxfixppsfpoind)))

    indxpsfp = arange(numbpsfp)
    indxpsfpinit = numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)

    dicttemp['indxfixpbacp'] = []
    for i in gdat.indxener:
        for c in indxback:
            dicttemp['indxfixpbacpene%dbac%d' % (i, c)] = cntr.incr()
            dicttemp['indxfixpbacp'].append(dicttemp['indxfixpbacpene%dbac%d' % (i, c)])
    dicttemp['indxfixpbacp'] = array(dicttemp['indxfixpbacp'])
    dicttemp['indxfixplenp'] = []
    dicttemp['indxfixpanglsour'] = []
    dicttemp['indxfixpanglhost'] = []
    dicttemp['indxfixpangllens'] = []
    dicttemp['indxfixpsour'] = []
    dicttemp['indxfixphost'] = []

    if gdat.pntstype == 'lens':
        dicttemp['indxfixplgalsour'] = cntr.incr()
        dicttemp['indxfixpbgalsour'] = cntr.incr()
        dicttemp['indxfixpspecsour'] = cntr.incr()
        dicttemp['indxfixpsizesour'] = cntr.incr()
        dicttemp['indxfixpellpsour'] = cntr.incr()
        dicttemp['indxfixpanglsour'] = cntr.incr()
        dicttemp['indxfixplgalhost'] = cntr.incr()
        dicttemp['indxfixpbgalhost'] = cntr.incr()
        dicttemp['indxfixpspechost'] = cntr.incr()
        dicttemp['indxfixpsizehost'] = cntr.incr()
        dicttemp['indxfixpbeinhost'] = cntr.incr()
        dicttemp['indxfixpellphost'] = cntr.incr()
        dicttemp['indxfixpanglhost'] = cntr.incr()
        dicttemp['indxfixpsherhost'] = cntr.incr()
        dicttemp['indxfixpsanghost'] = cntr.incr()
        dicttemp['indxfixpsour'] = [dicttemp['indxfixplgalsour'], dicttemp['indxfixpbgalsour'], dicttemp['indxfixpspecsour'], dicttemp['indxfixpsizesour'], \
						dicttemp['indxfixpellpsour'], dicttemp['indxfixpanglsour']]
        dicttemp['indxfixphost'] = [dicttemp['indxfixplgalhost'], dicttemp['indxfixpbgalhost'], dicttemp['indxfixpspechost'], dicttemp['indxfixpsizehost'], \
						dicttemp['indxfixpbeinhost'], dicttemp['indxfixpellphost'], dicttemp['indxfixpanglhost'], dicttemp['indxfixpsherhost'], dicttemp['indxfixpsanghost']]
        dicttemp['indxfixpemishost'] = [dicttemp['indxfixplgalhost'], dicttemp['indxfixpbgalhost'], dicttemp['indxfixpspechost'], dicttemp['indxfixpellphost'], \
						dicttemp['indxfixpanglhost']]
        dicttemp['indxfixplenp'] = list(set(dicttemp['indxfixpsour'] + dicttemp['indxfixphost'] + dicttemp['indxfixpemishost']))
    

    # number of fixed-dimension parameters
    numbfixp = cntr.incr(0)
    # indices of fixed-dimension parameters
    indxfixp = arange(numbfixp)

    ## number of model spectral parameters for each population
    numbspep = empty(numbpopl, dtype=int)
    liststrgspep = [[] for l in range(numbpopl)]
    liststrgfluxspep = [[] for l in range(numbpopl)]
    for l in range(numbpopl):
        if gdat.numbener > 1:
            liststrgspep[l] += ['sind']
            if spectype[l] == 'expo':
                liststrgspep[l] += ['expo']
            if spectype[l] == 'curv':
                liststrgspep[l] = ['curv']
        liststrgfluxspep[l] = ['flux'] + liststrgspep[l]
        numbspep[l] = len(liststrgspep[l]) 

    # number of components
    numbcomp = 3 + zeros(numbpopl, dtype=int)
    if gdat.numbener > 1:
        numbcomp += numbspep
    maxmnumbcomp = amax(numbcomp)

    indxcomp = []
    for l in indxpopl:
        indxcomp.append(arange(numbcomp[l]))

    # number of transdimensional parameters
    numbtrappopl = gdat.maxmnumbpnts * numbcomp
    numbtrapcumr = cumsum(numbtrappopl)
    numbtrapcuml = concatenate((array([0]), numbtrapcumr[:-1]))
    numbtrap = sum(numbtrappopl)
    
    # total number of parameters
    numbpara = numbfixp + numbtrap
    indxsampcomp = arange(numbfixp, numbpara)
    indxpara = arange(numbpara)

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
    
    for strg, valu in dicttemp.iteritems():
        
        if not isscalar(valu):
            continue
        k = valu
        
        namefixp[k] = strg[8:]

        if strg[:-1].endswith('pop'):
            
            if numbpopl == 1:
                strgpopl = ''
                strgpoplcomm = ''
            else:
                strgpopl = '%s' % strg[-1]
                strgpoplcomm = ',%s' % strg[-1]
            
            if namefixp[k].startswith('numbpnts'):
                strgfixp[k] = '$N_{%s}$' % strgpopl
                scalfixp[k] = 'pois'
                
            if namefixp[k].startswith('meanpnts'):
                strgfixp[k] = r'$\mu_{%s}$' % strgpopl
                scalfixp[k] = 'logt'
    
            if namefixp[k].startswith('gangdistscalp'):
                strgfixp[k] = r'$\gamma_{\theta%s}$' % strgpoplcomm
                scalfixp[k] = 'self'
            
            if namefixp[k].startswith('bgaldistscalp'):
                strgfixp[k] = r'$\gamma_{b%s}$' % strgpoplcomm
                scalfixp[k] = 'self'
            
            if namefixp[k].startswith('fluxdistslopp'):
                strgfixp[k] = r'$\alpha_{%s}$' % strgpopl
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('fluxdistbrek'):
                strgfixp[k] = '$f_{b%s}$' % strgpoplcomm
                scalfixp[k] = 'logt'
    
            if namefixp[k].startswith('fluxdistsloplowr'):
                strgfixp[k] = r'$\alpha_{l%s}$' % strgpoplcomm
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('fluxdistslopuppr'):
                strgfixp[k] = r'$\alpha_{u%s}$' % strgpoplcomm
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('sinddistmean'):
                strgfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablsind)
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('sinddiststdv'):
                strgfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablsind)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('curvdistmean'):
                strgfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablcurv)
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('curvdiststdv'):
                strgfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablcurv)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('expodistmean'):
                strgfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablexpo)
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('expodiststdv'):
                strgfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablexpo)
                scalfixp[k] = 'logt'
            
        if strg[:-1].endswith('evt'):
            
            if gdat.psfninfoprio:
                scalfixp[k] = 'gaus'
                n = k - dicttemp['indxfixpsigcene0evt0']
                meanfixp[k] = gdat.meanpsfp[n]
                stdvfixp[k] = gdat.meanpsfp[n]
            else:
                if strg[8:].startswith('sig'):
                    scalfixp[k] = 'logt'
                if strg[8:].startswith('gam'):
                    scalfixp[k] = 'atan'
                if strg[8:].startswith('psff'):
                    scalfixp[k] = 'atan'
                if strg[8:].startswith('onor'):
                    scalfixp[k] = 'logt'
                if strg[8:].startswith('oind'):
                    scalfixp[k] = 'atan'
                
            # strings for PSF parameters
            if strg[8:].startswith('sig'):
                strgvarbtemp = '\sigma'
                factfixpplot[k] = gdat.anglfact
            if strg[8:].startswith('gam'):
                strgvarbtemp = '\gamma'
            if strg[8:].startswith('psff'):
                strgvarbtemp = 'f'
            if strg[8:].startswith('onor'):
                strgvarbtemp = 'a'
            if strg[8:].startswith('oind'):
                strgvarbtemp = 'b'
            if strg[8:].startswith('sig') and psfntype == 'doubgaus' or psfntype == 'gausking' or psfntype == 'doubking':
                strgcomptemp = 'c'
            elif strg[8:].startswith('gam') and psfntype == 'gausking' or psfntype == 'doubking':
                strgcomptemp = 't'
            else:
                strgcomptemp = ''
            if gdat.numbener > 1:
                indxenertemp = gdat.indxenerincl[((k - dicttemp['indxfixpsigcene0evt0']) % (gdat.numbener * numbpsfptotl)) // numbpsfptotl]
                strgenertemp = '%s' % indxenertemp
            else:
                strgenertemp = ''
            if gdat.numbevtt > 1:
                indxevtttemp = gdat.indxevttincl[(k - dicttemp['indxfixpsigcene0evt0']) // (gdat.numbener * numbpsfptotl)]
                strgevtttemp = '%s' % indxevtttemp
            else:
                strgevtttemp = ''
            strgfixp[k] = r'$%s^{%s}_{%s%s}$' % (strgvarbtemp, strgcomptemp, strgenertemp, strgevtttemp)
        
        if strg[8:].startswith('bacp'):
            c = (k - dicttemp['indxfixpbacpene0bac0']) // gdat.numbener
            if gdat.numbener > 1:
                i = (k - dicttemp['indxfixpbacpene0bac0']) % gdat.numbener
                strgenertemp = '%d' % i
            else:
                strgenertemp = ''

            if numbback > 1:
                strgbacktemp = '%d' % c
            else:
                strgbacktemp = ''
            namefixp[k] = 'bacpene%dbac%d' % (i, c)
            strgfixp[k] = '$A_{%s%s}$' % (strgbacktemp, strgenertemp)
            scalfixp[k] = 'logt'
        
        if gdat.pntstype == 'lens':
            if k in dicttemp['indxfixplenp']:
                if strg[8:] == 'lgalsour':
                    strgfixp[k] = '$l_s$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strg[8:] == 'bgalsour':
                    strgfixp[k] = '$b_s$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strg[8:] == 'specsour':
                    strgfixp[k] = '$f_s$'
                    scalfixp[k] = 'logt'
                if strg[8:] == 'sizesour':
                    strgfixp[k] = '$a_s$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                if strg[8:] == 'ellpsour':
                    strgfixp[k] = r'$\epsilon_s$'
                    scalfixp[k] = 'self'
                if strg[8:] == 'anglsour':
                    strgfixp[k] = r'$\phi_s$'
                    scalfixp[k] = 'self'
                if strg[8:] == 'lgalhost':
                    strgfixp[k] = '$l_h$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strg[8:] == 'bgalhost':
                    strgfixp[k] = '$b_h$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strg[8:] == 'spechost':
                    strgfixp[k] = '$f_h$'
                    scalfixp[k] = 'logt'
                if strg[8:] == 'sizehost':
                    strgfixp[k] = '$a_h$'
                    scalfixp[k] = 'logt'
                if strg[8:] == 'beinhost':
                    strgfixp[k] = r'$\theta_{E,h}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                if strg[8:] == 'ellphost':
                    strgfixp[k] = r'$\epsilon_h$'
                    scalfixp[k] = 'self'
                if strg[8:] == 'anglhost':
                    strgfixp[k] = r'$\phi_h$'
                    scalfixp[k] = 'self'
                if strg[8:] == 'sherhost':
                    strgfixp[k] = r'$\gamma_e$'
                    scalfixp[k] = 'self'
                if strg[8:] == 'sanghost':
                    strgfixp[k] = r'$\phi_{\gamma}$'
                    scalfixp[k] = 'self'
        
        if scalfixp[k] == 'pois' or scalfixp[k] == 'self' or scalfixp[k] == 'logt' or scalfixp[k] == 'atan':
            
            if namefixp[k].startswith('numbpnts'):
                l = k
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k][:-4])[l]
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k][:-4])[l]
            elif namefixp[k][:-1].endswith('pop'):
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k][:-4])
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k][:-4])
            elif namefixp[k][:-1].endswith('evt') or namefixp[k][:-1].endswith('bac'):
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k][:-8])
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k][:-8])
            else:
                minmfixp[k] = getattr(gdat, 'minm' + namefixp[k])
                maxmfixp[k] = getattr(gdat, 'maxm' + namefixp[k])
        
        if scalfixp[k] == 'gaus' or scalfixp[k] == 'eerr':
            if gdat.psfninfoprio:
                meanfixp[k] = getattr(gdat, 'meanpsfp')[k-dicttemp['indxfixpsigcene0evt0']]
                stdvfixp[k] = getattr(gdat, 'stdvpsfp')[k-dicttemp['indxfixpsigcene0evt0']]
            else:
                if namefixp[k][:-1].endswith('pop'):
                    meanfixp[k] = getattr(gdat, 'mean' + namefixp[k][:-4])
                    stdvfixp[k] = getattr(gdat, 'stdv' + namefixp[k][:-4])
                elif namefixp[k][:-1].endswith('evt') or namefixp[k][:-1].endswith('bac'):
                    meanfixp[k] = getattr(gdat, 'mean' + namefixp[k][:-8])
                    stdvfixp[k] = getattr(gdat, 'stdv' + namefixp[k][:-8])
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
        
        if strg[16:20] == 'brek':
            strgfixpunit[k] = strgfixp[k] + ' [%s]' % getattr(gdat, 'lablfluxunit')
        elif strg[8:].startswith('sig'):
            strgfixpunit[k] = strgfixp[k] + ' [%s]' % getattr(gdat, 'lablgangunit')
        else:
            strgfixpunit[k] = strgfixp[k]

    namepara = list(deepcopy(namefixp))
    for k in range(numbtrap):
        indxpopltemp = (k - indxsampcomp[0]) // numbtrapcumr
        indxcomptemp = (k - indxsampcomp[0] - numbtrapcuml[indxpopltemp]) % numbcomp[indxpopl]
        namepara.append(liststrgcomp[indxpopltemp][indxcomptemp])
    
    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr':
            setattr(gdat, strgpara + strg, valu)

    for strg, valu in dicttemp.iteritems():
        # temp
        if isinstance(valu, ndarray) and strg[12:16] != 'dist':
            valu = sort(valu)
        setattr(gdat, strgpara + strg, valu)


def defn_truedefa(gdat, valu, strgvarb):
    
    strgtemp = 'trueindxfixp' + strgvarb
    for strg in gdat.__dict__.keys():
        if strg.startswith(strgtemp):
            strgrepl = 'true' + strgvarb + strg[len(strgtemp):]
            try:
                if getattr(gdat, strgrepl) == None:
                    setattr(gdat, strgrepl, valu)
            except:
                setattr(gdat, strgrepl, valu)


def setp_varbfull(gdat, strgpara, listfeat, typelimt='minmmaxm', numbpopl=None):
    
    numbfeat = len(listfeat)

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


def init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=None, indxevttplot=None, indxpoplplot=-1):

    if strg == 'this' or gdatmodi != None:
        pathfold = gdat.pathfram
    elif strg == 'true' or strg == '':
        pathfold = gdat.pathinit
    elif strg == 'post' or strg == 'medi' or strg == 'errr':
        pathfold = gdat.pathpost
    
    figr, axis = plt.subplots(figsize=(gdat.sizeimag, gdat.sizeimag))
    
    if indxenerplot == None:
        strgener = 'N'
    else:
        strgener = '%d' % gdat.indxenerincl[indxenerplot]
    
    if indxevttplot == None:
        strgevtt = 'N'
    elif indxevttplot == -1:
        strgevtt = 'A'
    else:
        strgevtt = '%d' % gdat.indxevttincl[indxevttplot]
    
    if indxpoplplot == -1:
        strgpopl = 'A'
    else:
        strgpopl = '%d' % indxpoplplot

    if gdatmodi == None:
        strgswep = ''
    else:
        strgswep = '_swep%09d' % gdatmodi.cntrswep
    
    path = '%s%s%s%s%s%s.pdf' % (pathfold, strgplot, strgener, strgevtt, strgpopl, strgswep)
   
    axis.set_xlabel(gdat.lablfeattotl['lgal'])
    axis.set_ylabel(gdat.lablfeattotl['bgal'])
    if indxenerplot != None and gdat.numbener > 1 or indxevttplot != None and gdat.numbevtt > 1:
        if indxenerplot != None and gdat.numbener > 1:
            titl = gdat.strgener[indxenerplot]
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


def retr_imag(gdat, axis, maps, strg, thisindxener=None, thisindxevtt=-1, cmap='Reds', vmin=None, vmax=None, scal=None, tdim=False):
    
    if scal == None:
        scal = gdat.scalmaps

    if vmin == None and vmax != None:
        vmin = -vmax
    
    draw_frambndr(gdat, axis)
    
    # flatten the array
    if tdim:
        if thisindxener == None:
            if strg == 'post':
                maps = maps.reshape((gdat.numbpixl, 4))
            else:
                maps = maps.reshape((gdat.numbpixl))
        else:
            if strg == 'post':
                maps = maps.reshape((gdat.numbener, gdat.numbpixl, gdat.numbevtt, 4))
            else:
                maps = maps.reshape((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    # take the relevant energy and PSF bins
    if thisindxener != None:
        if thisindxevtt == -1:
            maps = sum(maps[thisindxener, ...], axis=1)
        else:
            if strg == 'post':
                maps = maps[thisindxener, :, thisindxevtt, :]
            else:
                maps = maps[thisindxener, :, thisindxevtt]
    
    # project the map to 2D
    if gdat.pixltype == 'heal':
        maps = tdpy.util.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                            minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                            minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata)
    
    if gdat.pixltype == 'cart':
        shap = [gdat.numbsidecart] + list(maps.shape)
        shap[1] = gdat.numbsidecart
        maps = maps.reshape(shap).swapaxes(0, 1)
   
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

    axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, label='Model PS', marker='+', linewidth=2, color='b')
    
    if gdat.trueinfo:
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablhits, marker='x', linewidth=2, color='g')
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablbias, marker='*', linewidth=2, color='g', facecolor='none')
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, facecolor='none', \
                                                                                                label=gdat.truelablmiss, marker='o', linewidth=2, color='g')
    if gdat.pntstype == 'lens':
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='Model Source', marker='<', linewidth=2, color='b')

        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='Model Host', marker='s', linewidth=2, color='b')
        if gdat.trueinfo:
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='%s Source' % gdat.truelabl, marker='>', linewidth=2, color='g')
        
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='%s Host' % gdat.truelabl, marker='D', linewidth=2, color='g')
        
    axis.legend(bbox_to_anchor=[0.5, 1.1], loc='center', ncol=4)
        

def supr_fram(gdat, gdatmodi, axis, indxpoplplot=-1, trueonly=False):

    # true catalog
    if gdat.trueinfo:
       
        if indxpoplplot == -1:
            listindxpoplplot = gdat.indxpopl
        else:
            listindxpoplplot = [indxpoplplot]
        
        for l in listindxpoplplot:
            ## get the true catalog
            if gdat.numbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist, :].flatten())
                lgal = copy(gdat.truelgal[l])
                bgal = copy(gdat.truebgal[l])
                numbpnts = int(gdat.truenumbpnts[l])
                
                if gdatmodi != None and not trueonly:
                    
                    ## associations
                    ### missed
                    indx = gdatmodi.trueindxpntsassc[l].miss
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablmiss, facecolor='none', \
                                                                                                                                marker='o', linewidth=2, color='g')
                    
                    ### biased
                    indx = gdatmodi.trueindxpntsassc[l].bias
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablbias, marker='*', linewidth=2, color='g', facecolor='none')
                    
                    ### hit
                    indx = gdatmodi.trueindxpntsassc[l].hits
                    
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
                    mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]])
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
                                        gdat.fluxfactplot * gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l][k]], edgecolor='b', facecolor='none', ls='--', lw=2))
    
    
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
        indxpsfponor = numbpsfpform + numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)
        indxpsfpoind = numbpsfpform + numbpsfptotl * arange(gdat.numbener * gdat.numbevtt) + 1
    else:
        indxpsfponor = []
        indxpsfpoind = []

    return numbpsfpform, numbpsfpoaxi, numbpsfptotl, indxpsfponor, indxpsfpoind


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


def pert_llik(gdat, gdatmodi, indxparapert, stdvparapert):

    numbpert = indxparapert.size 
    gdatmodi.nextsamp = copy(gdatmodi.thissamp)
    for k in range(numbpert):
        gdatmodi.nextsamp[indxparapert[k]] += stdvparapert[k]
    
    gdatmodi.nextsampvarb = retr_sampvarb(gdat, gdatmodi.thisindxpntsfull, gdatmodi.nextsamp, 'next')
    
    lpospert = retr_negalpos(gdat, gdatmodi)
    
    return lpospert


def retr_deflextr(gdat, sher, sang):
    
    factcosi = sher * cos(2. * sang)
    factsine = sher * cos(2. * sang)
    defllgal = factcosi * gdat.lgalgridcart + factsine * gdat.bgalgridcart
    deflbgal = factsine * gdat.lgalgridcart - factcosi * gdat.bgalgridcart
    
    return dstack((defllgal, deflbgal)) 


def retr_defl(gdat, lgal, bgal, bein, ellp, angl, rcor, thisevalcirc=False):
    
    if thisevalcirc and gdat.evalcirc:
        indxfluxproxtemp = digitize(bein[k], gdat.binsfluxprox) - 1
        indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
        indxpixlprox = gdat.indxsidemeshprox[indxfluxproxtemp][indxpixlpnts]
    else:
        indxpixlprox = gdat.indxsidemesh
    
    # translate the grid
    lgaltran = gdat.lgalgridcart[indxpixlprox] - lgal
    bgaltran = gdat.bgalgridcart[indxpixlprox] - bgal
    
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
    defllgal = cos(angl) * defllgalrttr + sin(angl) * deflbgalrttr
    deflbgal = -sin(angl) * defllgalrttr + cos(angl) * deflbgalrttr
    
    
    return dstack((defllgal, deflbgal))
   

def retr_negalpos(gdat, gdatmodi):
   
    proc_samp(gdat, gdatmodi, 'next')
    
    return -gdatmodi.nextlliktotl - gdatmodi.nextlpritotl


def retr_lprifluxdist(gdat, gdatmodi, flux, sampvarb, l):
    
    if gdat.fluxdisttype[l] == 'powr':
        fluxdistslop = sampvarb[gdat.indxfixpfluxdistslop[l]]
        lprifluxdist = sum(log(pdfn_flux_powr(gdat, flux, fluxdistslop)))
    if gdat.fluxdisttype[l] == 'brok':
        fluxdistbrek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
        fluxdistsloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
        fluxdistslopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
        lprifluxdist = sum(log(pdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)))

    return lprifluxdist


def retr_lprisinddist(gdat, gdatmodi, sind, sampvarb, l):
    
    lprisinddist = sum(log(pdfn_gaus(sind, sampvarb[gdat.indxfixpsinddistmean[l]], sampvarb[gdat.indxfixpsinddiststdv[l]]))) 
    
    return lprisinddist


def retr_lpricurvdist(gdat, gdatmodi, curv, sampvarb, l):
    
    lpri = sum(log(pdfn_gaus(curv, sampvarb[gdat.indxfixpcurvdistmean[l]], sampvarb[gdat.indxfixpcurvdiststdv[l]])))
    
    return lpri


def retr_lpriexpodist(gdat, gdatmodi, expo, sampvarb, l):
    
    lpri = sum(log(pdfn_gaus(expo, sampvarb[gdat.indxfixpexpodistmean[l]], sampvarb[gdat.indxfixpexpodiststdv[l]])))
    
    return lpri


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
        fluxdisttype = gdat.truefluxdisttype
        spatdisttype = gdat.truespatdisttype
    else:
        strgtype = ''
        spectype = gdat.spectype
        fluxdisttype = gdat.fluxdisttype
        spatdisttype = gdat.spatdisttype
    
    # common dictionary
    dicttemp = {}
           
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
        defl = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost, 0.)
        deflextr = retr_deflextr(gdat, sherhost, sanghost)
        defl += deflextr

    if gdat.pntstype == 'lght':
        varioaxi = getattr(gdat, strgtype + 'varioaxi')

    indxpntsfull = list(getattr(gdatobjt, strg + 'indxpntsfull'))
    indxsamplgal, indxsampbgal, indxsampflux, indxsampsind, indxsampcurv, indxsampexpo, indxsampcomp = retr_indx(gdat, indxpntsfull, spectype)
    dicttemp['indxsamplgal'], dicttemp['indxsampbgal'], dicttemp['indxsampflux'], dicttemp['indxsampsind'], dicttemp['indxsampcurv'], \
																				dicttemp['indxsampexpo'], dicttemp['indxsampcomp'] = retr_indx(gdat, indxpntsfull, spectype)
    for strgcomp in gdat.liststrgcompdefa:
        setattr(gdatobjt, strg + 'indxsamp' + strgcomp, dicttemp['indxsamp' + strgcomp])
    setattr(gdatobjt, strg + 'indxsampcomp', dicttemp['indxsampcomp'])
    
    if strg == 'next' and gdat.verbtype > 1:
        show_samp(gdat, gdatmodi)
    
    numbpnts = getattr(gdatobjt, strg + 'sampvarb')[gdat.indxfixpnumbpnts].astype(int)
    numbpopl = numbpnts.size
    
    for strgcomp in gdat.liststrgcompdefa:
        dicttemp[strgcomp] = [[] for l in range(numbpopl)] 
    for l in range(numbpopl):
    	for strgcomp in gdat.liststrgcomp[l]:
            dicttemp[strgcomp][l] = sampvarb[dicttemp['indxsamp' + strgcomp][l]]
    
    dicttemp['spec'] = [[] for l in range(numbpopl)]
    for l in range(numbpopl):
        dicttemp['spec'][l] = retr_spec(gdat, dicttemp['flux'][l], dicttemp['sind'][l], dicttemp['curv'][l], dicttemp['expo'][l], spectype=spectype[l])
        
    lgalconc = concatenate(dicttemp['lgal'])
    bgalconc = concatenate(dicttemp['bgal'])
    specconc = concatenate(dicttemp['spec'], axis=1)
    numbpntsconc = lgalconc.size
    
    # process a sample vector and the occupancy list to calculate secondary variables
	## secondary variables needed for likelihood evaluation    
    psfp = sampvarb[getattr(gdat, 'indxfixppsfp')]
    if gdat.pntstype == 'lens':
        psfnkern = AiryDisk2DKernel(psfp[0] / gdat.sizepixl)
        
    if gdat.pntstype == 'lght':
        ### PSF off-axis factor
        if varioaxi:
            onor = sampvarb[getattr(gdat, 'indxfixppsfponor')]
            oind = sampvarb[getattr(gdat, 'indxfixppsfpoind')]
            factoaxi = retr_factoaxi(gdat, gdat.binsoaxi, onor, oind)
    
        psfntype = getattr(gdat, strgtype + 'psfntype')
        psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, varioaxi)

        if varioaxi:
            psfnintp = []
            for p in gdat.indxoaxi:
                psfnintp.append(interp1d_pick(gdat.binsangl, psfn[:, :, :, p], axis=1))
        else:
            psfnintp = interp1d_pick(gdat.binsangl, psfn, axis=1)
        setattr(gdatobjt, strg + 'psfnintp', psfnintp)
   
    bacp = sampvarb[getattr(gdat, 'indxfixpbacp')]
    
    if gdat.pntstype == 'lens':
        
        ## components
        if numbpntsconc > 0 and not raww:
            for k in range(numbpntsconc):
                # temp -- fix truncation
                defl += retr_defl(gdat, lgalconc[k], bgalconc[k], specconc[0, k], 0., 0., 0., thisevalcirc=True)
                        
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
            lpri = zeros(gdat.numblpri)
            for l in gdat.indxpopl:
                lpri[0] -= 0.5 * gdat.priofactdoff * gdat.numbcomp[l] * numbpnts[l]
                lpri[1+0*gdat.numbpopl+l] = retr_probpois(numbpnts[l], meanpnts[l])
                if spatdisttype[l] == 'unif':
                    lpri[1+gdat.numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang)
                    lpri[1+2*gdat.numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang)
                elif spatdisttype[l] == 'disc':
                    lpri[1+gdat.numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang)
                    lpri[1+2*gdat.numbpopl+l] = sum(log(pdfn_dexp(dicttemp['bgal'][l], gdat.maxmgang, sampvarb[gdat.indxfixpbgaldistscal[l]]))) 
                    if not isfinite(sum(log(pdfn_dexp(dicttemp['bgal'][l], gdat.maxmgang, sampvarb[gdat.indxfixpbgaldistscal[l]])))):
                        print 'dicttemp[bgal][l]'
                        print dicttemp['bgal'][l]
                elif spatdisttype[l] == 'gang':
                    gang = retr_gang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                    lpri[1+gdat.numbpopl+l] = sum(log(pdfn_expo(gang, gdat.maxmgang, sampvarb[gdat.indxfixpgangdistscal[l]]))) 
                    #if not isfinite(sum(log(pdfn_expo(gang, gdat.maxmgang, sampvarb[gdat.indxfixpgangdistscal[l]])))):  
                    #    print 'gang'
                    #    print gang

                    # temp
                    #lpri[1+2*gdat.numbpopl+l] = -numbpnts[l] * log(2. * pi) 
                    lpri[1+2*gdat.numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang) 
                
                lpri[1+3*gdat.numbpopl+l] = retr_lprifluxdist(gdat, gdatmodi, dicttemp['flux'][l], sampvarb, l)
                if gdat.numbener > 1:
                    lpri[1+4*gdat.numbpopl+l] = retr_lprisinddist(gdat, gdatmodi, dicttemp['sind'][l], sampvarb, l)
                    if gdat.spectype[l] == 'curv':
                        lpri[1+5*gdat.numbpopl+l] = retr_lpricurvdist(gdat, gdatmodi, dicttemp['curv'][l], sampvarb, l)
                    if gdat.spectype[l] == 'expo':
                        lpri[1+6*gdat.numbpopl+l] = retr_lpriexpodist(gdat, gdatmodi, dicttemp['expo'][l], sampvarb, l)
            
            lpritotl = sum(lpri)
            setattr(gdatmodi, strg + 'lpri', lpri)
            
            if strg == 'next' and (gdatmodi.propbrth or gdatmodi.propdeth):
                
                if gdatmodi.propbrth:
                    sampvarbtemp = gdatmodi.nextsampvarb
                if gdatmodi.propdeth:
                    sampvarbtemp = gdatmodi.thissampvarb
                
                gdatmodi.thislpau = zeros(gdat.numblpau)
                gdatmodi.thislpau[gdat.maxmnumbcomp*l+0] = -log(2. * gdat.maxmgang)
                gdatmodi.thislpau[gdat.maxmnumbcomp*l+1] = -log(2. * gdat.maxmgang)
                gdatmodi.thislpau[gdat.maxmnumbcomp*l+2] = retr_lprifluxdist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[2]], gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                if gdat.numbener > 1:
                    gdatmodi.thislpau[gdat.maxmnumbcomp*l+3] = retr_lprisinddist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[3]], \
                                                                                                                gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'curv':
                        gdatmodi.thislpau[gdat.maxmnumbcomp*l+4] = retr_lpricurvdist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[4]], \
                                                                                                                gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'expo':
                        gdatmodi.thislpau[gdat.maxmnumbcomp*l+4] = retr_lpriexpodist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[4]], \
                                                                                                                gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)

                if gdatmodi.propbrth:
                    gdatmodi.thislpau *= -1.
                
                gdatmodi.thislpautotl = sum(gdatmodi.thislpau)
                
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
            setattr(gdatmodi, 'thisdeltllik', deltllik)
       
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
        setattr(gdatmodi, strg + 'llik', gdatmodi.templlik) 
    
    if strg != 'next':
        
        # temp
        #setattr(gdatobjt, strg + 'spec', dicttemp['spec'])
        
        if strg == 'this':
            lpostotl = lpritotl + gdatmodi.templliktotl
            setattr(gdatobjt, strg + 'lpostotl', lpostotl) 
       
        setattr(gdatobjt, strg + 'psfp', psfp)
        if gdat.pixltype != 'unbd':
            resicnts = gdat.datacnts - modlcnts
            setattr(gdatobjt, strg + 'resicnts', resicnts)
            
        ## component features
        if gdat.numbtrap > 0:
            
            ### derived quantities
            dicttemp['gang'] = []
            dicttemp['aang'] = []
            dicttemp['cnts'] = []
            for l in range(numbpopl):
                #### radial and angular coordinates
                dicttemp['gang'].append(retr_gang(dicttemp['lgal'][l], dicttemp['bgal'][l]))
                dicttemp['aang'].append(retr_aang(dicttemp['lgal'][l], dicttemp['bgal'][l]))
                
                #### number of expected counts
                dicttemp['cnts'].append(retr_pntscnts(gdat, dicttemp['lgal'][l], dicttemp['bgal'][l], dicttemp['spec'][l]))

            ### distribution of component features and derived quantities
            if gdat.numbener > 1:
                dicttemp['fluxsindhist'] = zeros((numbpopl, gdat.numbbinsplot, gdat.numbbinsplot))
            for strgfeat in gdat.liststrgfeat:
                if strgfeat == 'spec':
                    temp = zeros((numbpopl, gdat.numbbinsplot, gdat.numbener))
                elif strgfeat == 'cnts':
                    temp = zeros((numbpopl, gdat.numbbinsplot, gdat.numbener, gdat.numbevtt))
                else:
                    temp = zeros((numbpopl, gdat.numbbinsplot))
                dicttemp[strgfeat + 'hist'] = temp
           
            ## feature distributions of elements
            # temp -- this should be done for true sources
            # find the indices of the model PSs that are in the comparison area
            indxmodlpntscomp = [[] for l in range(numbpopl)]
            for l in range(numbpopl):
                indxmodlpntscomp[l] = where((fabs(dicttemp['lgal'][l]) < gdat.maxmgangdata) & (fabs(dicttemp['bgal'][l]) < gdat.maxmgangdata))[0]
            
            ## histograms of element features
            ### flux - color
            if gdat.numbener > 1:
                for l in range(numbpopl):
                    dicttemp['fluxsindhist'][l, :, :] = histogram2d(dicttemp['flux'][l][indxmodlpntscomp[l]], dicttemp['sind'][l][indxmodlpntscomp[l]], \
                                                                                                                                [gdat.binsfluxplot, gdat.binssindplot])[0]
                setattr(gdatobjt, strg + 'fluxsindhist', dicttemp['fluxsindhist'])
            
            for strgfeat in gdat.liststrgfeat:
                for l in range(numbpopl):
                    if strgfeat == 'spec':
                        for i in gdat.indxener:
                            dicttemp[strgfeat + 'hist'][l, :, i] = histogram(dicttemp['spec'][l][i, indxmodlpntscomp[l]], gdat.binsspecplot[i, :])[0]
                    elif strgfeat == 'cnts':
                        for i in gdat.indxener:
                            for m in gdat.indxevtt:
                                dicttemp[strgfeat + 'hist'][l, :, i, m] = histogram(dicttemp['cnts'][l][i, indxmodlpntscomp[l], m], gdat.binscnts[i, :])[0]
                    elif not (strgfeat == 'curv' and spectype[l] != 'curv' or strgfeat == 'expo' and spectype[l] != 'expo'):
                        bins = getattr(gdat, 'bins' + strgfeat + 'plot')
                        dicttemp[strgfeat + 'hist'][l, :] = histogram(dicttemp[strgfeat][l][indxmodlpntscomp[l]], bins)[0]

                setattr(gdatobjt, strg + strgfeat + 'hist', dicttemp[strgfeat + 'hist'])
   
            ### priors on component features
            for strgfeat in gdat.liststrgfeat:
                dicttemp[strg + strgfeat + 'histprio'] = empty((numbpopl, gdat.numbbinsplotprio))
            
                # temp -- this does not work for mismodeling, need strg +
                minm = getattr(gdat, 'minm' + strgfeat + 'plot')
                maxm = getattr(gdat, 'maxm' + strgfeat + 'plot')
                for l in range(numbpopl):
                    if strgfeat in gdat.liststrgfeatprio[l]:
                        
                        xdat = getattr(gdat, 'mean' + strgfeat + 'plotprio')
                        deltprio = getattr(gdat, 'delt' + strgfeat + 'plotprio')
                        delt = getattr(gdat, 'delt' + strgfeat + 'plot')
                        
                        booltemp = False
                        if strgfeat == 'gang' and spatdisttype[l] == 'gang' or strgfeat == 'bgal' and spatdisttype[l] == 'disc': 
                            scal = sampvarb[getattr(gdat, 'indxfixp' + strgfeat + 'distscal')[l]]
                            if strgfeat == 'gang' and spatdisttype[l] == 'gang':
                                pdfn = pdfn_expo(xdat, maxm, scal)
                            else:
                                pdfn = pdfn_dexp(xdat, maxm, scal)
                            booltemp = True
                        if strgfeat == 'aang' and spatdisttype[l] == 'gang':
                            pdfn = 1. / 2. / pi + zeros_like(xdat)
                            booltemp = True
                        if spatdisttype[l] == 'unif' and (strgfeat == 'lgal' or strgfeat == 'bgal') or spatdisttype[l] == 'disc' and strgfeat == 'lgal': 
                            pdfn = 1. / 2. / gdat.maxmgang + zeros_like(xdat)
                            booltemp = True
                        if strgfeat == 'flux' and fluxdisttype[l] == 'powr':
                            slop = sampvarb[getattr(gdat, 'indxfixp' + strgfeat + 'distslop')[l]]
                            pdfn = pdfn_powr(xdat, slop, minm, maxm)
                            booltemp = True
                        elif strgfeat == 'flux' and fluxdisttype[l] == 'brok':
                            brek = sampvarb[gdat.indxfixpfluxdistbrek[l]]
                            sloplowr = sampvarb[gdat.indxfixpfluxdistsloplowr[l]]
                            slopuppr = sampvarb[gdat.indxfixpfluxdistslopuppr[l]]
                            pdfn = pdfn_brok(xdat, brek, sloplowr, slopuppr, minm, maxm)
                            booltemp = True
                        elif strgfeat == 'sind' or strgfeat == 'curv' and spectype[l] == 'curv' or strgfeat == 'expo' and spectype[l] == 'expo':
                            # this does not work for mismodeling
                            mean = sampvarb[getattr(gdat, 'indxfixp' + strgfeat + 'distmean')[l]]
                            stdv = sampvarb[getattr(gdat, 'indxfixp' + strgfeat + 'diststdv')[l]]
                            if strgfeat == 'expo' and spectype[l] == 'expo':
                                pdfn = pdfn_gaus(xdat, mean, stdv)
                            else:
                                pdfn = pdfn_gaus(xdat, mean, stdv)
                            booltemp = True
                        
                        if booltemp:
                            dicttemp[strg + strgfeat + 'histprio'][l, :] = meanpnts[l] * pdfn * deltprio * delt[0] / deltprio[0]
                
                setattr(gdatobjt, strg + strgfeat + 'histprio', dicttemp[strg + strgfeat + 'histprio'])
        
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
                sigm = []
                for l in gdat.indxpopl:
                    # temp -- zero exposure pixels will give zero counts
                    if gdat.varioaxi:
                        sigmtemp = retr_sigm(gdat, dicttemp['cnts'][l], cntsbackfwhm, lgal=lgal[l], bgal=bgal[l])
                    else:
                        sigmtemp = retr_sigm(gdat, dicttemp['cnts'][l], cntsbackfwhm)
                    sigm.append(sigmtemp)
                
            if gdat.calcerrr and gdat.numbtrap > 0:
                pntsflux = retr_pntsflux(gdat, lgalconc, bgalconc, specconc, psfnintp, gdat.varioaxi, evalcirc=False)
                pntscnts = retr_cntsmaps(gdat, pntsflux)
                errrcnts = pntscnts - temppntscnts
                indxcubegood = where(temppntscnts > 1e-10)
                setattr(gdatobjt, strg + 'errrcnts', errrcnts)
                if False and amax(fabs(errrcnts)) > 0.1:
                    raise Exception('Approximation error in calculating the PS flux map is above the tolerance level.')

            #fluxbrgt, fluxbrgtassc = retr_fluxbrgt(gdat, lgalconc, bgalconc, concatenate(dicttemp['flux']))
            #setattr(gdatobjt, strg + 'fluxbrgt', fluxbrgt)
            #setattr(gdatobjt, strg + 'fluxbrgtassc', fluxbrgtassc)

        if gdat.pntstype == 'lens':
            lenscnts = retr_cntsmaps(gdat, lensflux, cart=True)
            hostcntsmaps = retr_cntsmaps(gdat, hostfluxmaps, cart=True)
            
            setattr(gdatobjt, strg + 'psfnkern', psfnkern)
            setattr(gdatobjt, strg + 'modlfluxuncv', modlfluxuncv)
            setattr(gdatobjt, strg + 'lenscnts', lenscnts)
            setattr(gdatobjt, strg + 'hostcntsmaps', hostcntsmaps)

            ### deflection
            #### number of deflection components
            if numbpntsconc > 0:
                indxpntssortbrgt = argsort(dicttemp['flux'][0])[::-1]
                lgalsort = dicttemp['lgal'][0][indxpntssortbrgt][:gdat.numbdeflpnts]
                bgalsort = dicttemp['bgal'][0][indxpntssortbrgt][:gdat.numbdeflpnts]
                beinsort = dicttemp['flux'][0][indxpntssortbrgt][:gdat.numbdeflpnts]
            
            deflsing = zeros((gdat.numbsidecart, gdat.numbsidecart, 2, gdat.numbdeflsing))
            for k in range(gdat.numbdeflsing):
                if k == 0:
                    deflsing[:, :, :, k] = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost, 0.)
                elif k == 1:
                    deflsing[:, :, :, k] = deflextr
                else:
                    # temp
                    if k - 2 < lgalsort.size:
                        deflsing[:, :, :, k] = retr_defl(gdat, lgalsort[k-2], bgalsort[k-2], beinsort[k-2], 0., 0., 0.)

            ### convergence
            conv = retr_conv(gdat, defl) 
            convpsec = retr_psec(gdat, conv)
            convpsecodim = retr_psecodim(gdat, convpsec) 
            
            histdefl = histogram(defl, bins=gdat.binsdeflplot)[0]
            setattr(gdatobjt, strg + 'conv', conv)
            setattr(gdatobjt, strg + 'convpsec', convpsec)
            setattr(gdatobjt, strg + 'convpsecodim', convpsecodim)
            setattr(gdatobjt, strg + 'histdefl', histdefl)
            setattr(gdatobjt, strg + 'deflsing', deflsing)
     
	    ### PS indices to compare with the reference catalog
        if strg == 'this' and gdat.numbtrap > 0 and gdat.trueinfo:
            
            if gdat.pntstype == 'lens' and gdat.trueinfo and gdat.datatype == 'mock':
                gdatmodi.thisdeflsingresi = gdatmodi.thisdeflsing - gdat.truedeflsing
                gdatmodi.thisdeflresi = gdatmodi.thisdefl - gdat.truedefl
                gdatmodi.thisdeflcomp = 1. - sum(gdatmodi.thisdefl * gdat.truedefl, axis=2) / sqrt(sum(gdatmodi.thisdefl**2, axis=2)) / sqrt(sum(gdat.truedefl**2, axis=2))
           
        if gdatmodi != None:
            ### corrected prior
            gdatmodi.thislprinorm = 0.
            for l in gdat.indxpopl:
                # temp -- brok terms are not complete
                break
                numbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]]
                meanpnts = gdatmodi.thissampvarb[gdat.indxfixpmeanpnts[l]]
                gdatmodi.thislprinorm += numbpnts * gdat.priofactlgalbgal + gdat.priofactfluxdistslop + gdat.priofactmeanpnts - log(meanpnts)
                flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]]
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

    # correlate the current catalog sample with the reference catalog
    if strg == 'this':
        gdatmodi.trueindxpntsassc = [tdpy.util.gdatstrt() for l in gdat.trueindxpopl]
        gdatmodi.thisspecassc = [[] for l in gdat.trueindxpopl] 
        for l in gdat.trueindxpopl:
            gdatmodi.trueindxpntsassc[l].miss = []
            gdatmodi.trueindxpntsassc[l].bias = []
            gdatmodi.trueindxpntsassc[l].hits = []
            gdatmodi.trueindxpntsassc[l].mult = []
            if gdat.asscmetrtype == 'dist':
                dir1 = array([gdat.truelgal[l], gdat.truebgal[l]])
    
            indxmodlpnts = zeros(gdat.truenumbpnts[l], dtype=int) - 1
            specassc = zeros((gdat.numbener, gdat.truenumbpnts[l]))
            numbassc = zeros(gdat.truenumbpnts[l])
            metrassc = zeros(gdat.truenumbpnts[l]) + 3 * gdat.maxmgang
        
            for k in range(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]].astype(int)):
               
                # determine which true PSs satisfy the match criterion
                if gdat.asscmetrtype == 'dist':
                    dir2 = array([dicttemp['lgal'][l][k], dicttemp['bgal'][l][k]])
                    metr = retr_angldist(gdat, dir1, dir2)
                    trueindxpntstemp = where(metr < gdat.anglassc)[0]
                if gdat.asscmetrtype == 'defl':
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
                        specassc[:, trueindxpntstemp[0]] = dicttemp['spec'][l][:, k]
                        metrassc[trueindxpntstemp[0]] = metr[0]
                        indxmodlpnts[trueindxpntstemp[0]] = k
    
            # get the flux limit that delineates the biased associations and hits 
            fluxbias = empty((2, gdat.numbener, gdat.truenumbpnts[l]))
            for i in gdat.indxener:
                fluxbias[:, i, :] = retr_fluxbias(gdat, gdat.truespec[l][0, i, :], i)
    
            # divide associations into subgroups
            for k in range(gdat.truenumbpnts[l]):
                if numbassc[k] == 0:
                    gdatmodi.trueindxpntsassc[l].miss.append(k)
                else:
                    if numbassc[k] > 1:
                        gdatmodi.trueindxpntsassc[l].mult.append(k)
            
                    ## check whether the flux of the associated model point source matches well with the flux of the deterministic point source
                    for i in gdat.indxener:
                        boolbias = specassc[i, k] > fluxbias[1, i, k] or specassc[i, k] < fluxbias[0, i, k]
                        if boolbias:
                            gdatmodi.trueindxpntsassc[l].bias.append(k)
                        else:
                            gdatmodi.trueindxpntsassc[l].hits.append(k)
       
            gdatmodi.thisspecassc[l] = zeros((gdat.numbener, gdat.truenumbpnts[l]))
            temp = where(indxmodlpnts >= 0)[0]
            gdatmodi.thisspecassc[l][:, temp] = dicttemp['spec'][l][:, indxmodlpnts[temp]]
    
    
def retr_indxpntscomp(gdat, lgal, bgal):


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


