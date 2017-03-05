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
        oaxitype = True
        numboaxi = psfn.shape[3]
        wdth = zeros((gdat.numbener, gdat.numbevtt, numboaxi))
    else:
        oaxitype = False
        numboaxi = psfn.shape[2]
        wdth = zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            for p in arange(numboaxi):
                if oaxitype:
                    psfntemp = psfn[i, :, m, p]
                else:
                    if p > 0:
                        break
                    psfntemp = psfn[i, :, m]
                indxanglgood = argsort(psfntemp)
                intpwdth = max(frac * amax(psfntemp), amin(psfntemp))
                if intpwdth > amin(psfntemp[indxanglgood]) and intpwdth < amax(psfntemp[indxanglgood]):
                    wdthtemp = interp1d_pick(psfntemp[indxanglgood], gdat.binsanglplot[indxanglgood])(intpwdth)
                else:
                    wdthtemp = 0.
                if oaxitype:
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


def retr_indxsampcomp(gdat, indxpntsfull, spectype):
    
    indxsamplgal = [[] for l in gdat.indxpopl]
    indxsampbgal = [[] for l in gdat.indxpopl]
    indxsampflux = [[] for l in gdat.indxpopl]
    indxsampsind = [[] for l in gdat.indxpopl]
    indxsampcurv = [[] for l in gdat.indxpopl]
    indxsampexpo = [[] for l in gdat.indxpopl]
    indxsampcomp = []
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            indxsamptemp = gdat.indxsampcompinit + gdat.numbtrapcuml[l] + array(indxpntsfull[l], dtype=int) * gdat.numbcomp[l]
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


def retr_pntsflux(gdat, lgal, bgal, spec, psfnintp, oaxitype, evalcirc):
    
    if lgal.ndim == 0:
        lgal = array([lgal])
        bgal = array([bgal])
        spec = spec[:, None]

    numbpnts = lgal.size
    pntsfluxsing = zeros((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    for k in range(numbpnts):
        
        if evalcirc:
            indxfluxproxtemp = digitize(spec[gdat.indxenerfluxdist[0], k], gdat.binsprox) - 1
            indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
           
            # temp
            #if gdat.strgcnfg == 'pcat_ferm_quas_mock':
            #    print 'spec'
            #    print spec
            #    print spec.shape
            #    print 'spec[gdat.indxenerfluxdist[0], k]'
            #    print spec[gdat.indxenerfluxdist[0], k]
            #    print 'gdat.binsprox'
            #    print gdat.binsprox
            #    print 'digitize(spec[gdat.indxenerfluxdist[0], k], gdat.binsprox)'
            #    print digitize(spec[gdat.indxenerfluxdist[0], k], gdat.binsprox)
            #    print 'indxfluxproxtemp'
            #    print indxfluxproxtemp
            #    print 'indxpixlpnts'
            #    print indxpixlpnts
            #    summgene(indxpixlpnts)
            #    print

            indxpixltemp = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
            if isinstance(indxpixltemp, int):
                indxpixltemp = gdat.indxpixl
        else:
            indxpixltemp = gdat.indxpixl
    
        # calculate the distance to all pixels from each point source
        dist = retr_angldistunit(gdat, lgal[k], bgal[k], indxpixltemp)
        
        if oaxitype:
            indxoaxitemp = retr_indxoaxipnts(gdat, lgal[k, None], bgal[k, None])
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
        if gdat.specback[c] != None:
            norm = gdat.specback[c] * bacp[gdat.indxbacpback[c]]
        else:
            norm = bacp[gdat.indxbacpback[c]]
        modlflux += norm[:, None, None] * gdat.backflux[c][tempindx]        

    return modlflux


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


def pdfn_powr(xdat, minm, maxm, slop):
  
    norm = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop))
    
    pdfn = norm * xdat**(-slop)
    
    return pdfn


def pdfn_bind(xdat, minm, maxm, bins, norm):
  
    pdfn = norm / trapz(norm, bins)
    pdfn = interp(xdat, bins, pdfn)

    return pdfn


def cdfn_bind(para, minm, maxm, bins, norm):
  
    pdfn = norm / trapz(norm, bins)
    cdfn = cumsum(diff(bins) * (pdfn[:-1] + pdfn[1:]) / 2.)
    cdfn = concatenate((array([0.]), cdfn))
    paraunit = interp(para, bins, cdfn)

    return paraunit


def icdf_bind(paraunit, minm, maxm, bins, norm):
  
    pdfn = norm / trapz(norm, bins)
    cdfn = cumsum(diff(bins) * (pdfn[:-1] + pdfn[1:]) / 2.)
    cdfn = concatenate((array([0.]), cdfn))
    para = interp(paraunit, cdfn, bins)

    return para


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


def cdfn_atan(para, minmpara, maxmpara):
    
    paraunit = (arctan(para) - arctan(minmpara)) / (arctan(maxmpara) - arctan(minmpara))
    
    return paraunit


def icdf_atan(paraunit, minmpara, maxmpara):

    para = tan((arctan(maxmpara) - arctan(minmpara)) * paraunit + arctan(minmpara))
    
    return para


def pdfn_atan(para, minmpara, maxmpara):

    pdfn = 1. / (para**2 + 1.) / (arctan(maxmpara) - arctan(minmpara))
    
    return pdfn


def cdfn_fixp(gdat, strg, fixp, thisindxfixp):
    
    scalfixp = getattr(gdat, strg + 'scalfixp')[thisindxfixp]
    
    if scalfixp == 'self' or scalfixp == 'logt' or scalfixp == 'atan':
        
        minmfixp = getattr(gdat, strg + 'minmfixp')[thisindxfixp]
        factfixp = getattr(gdat, strg + 'factfixp')[thisindxfixp]

        if scalfixp == 'self':
            fixpunit = cdfn_self(fixp, minmfixp, factfixp)
        elif scalfixp == 'logt':
            fixpunit = cdfn_logt(fixp, minmfixp, factfixp)
        elif scalfixp == 'atan':
            maxmfixp = getattr(gdat, strg + 'maxmfixp')[thisindxfixp]
            fixpunit = cdfn_atan(fixp, minmfixp, maxmfixp)
    
    elif scalfixp == 'gaus' or scalfixp == 'eerr':
        meanfixp = getattr(gdat, strg + 'meanfixp')[thisindxfixp]
        stdvfixp = getattr(gdat, strg + 'stdvfixp')[thisindxfixp]
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
            maxmfixp = getattr(gdat, strg + 'maxmfixp')[thisindxfixp]
            fixp = icdf_atan(fixpunit, minmfixp, maxmfixp)
    
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


def retr_thisindxprop(gdat, gdatmodi, thisindxpopl=None, brth=False, deth=False):

    gdatmodi.thisaccppsfn = True
    gdatmodi.thisaccpprio = True
    
    # initialize the Boolean flag indicating the type of transdimensional proposal
    gdatmodi.propwith = False
    gdatmodi.propbrth = False
    gdatmodi.propdeth = False
    gdatmodi.propsplt = False
    gdatmodi.propmerg = False
    
    gdatmodi.thisjcbnfact = 0.
    gdatmodi.thislpautotl = 0. 
    gdatmodi.thiscombfact = 0.
   
    # index of the population in which a transdimensional proposal will be made
    if thisindxpopl == None:
        gdatmodi.indxpoplmodi = choice(gdat.indxpopl)
    else:
        gdatmodi.indxpoplmodi = thisindxpopl

    if (rand() < gdat.probtran or brth or deth) and (gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] != gdat.minmnumbpnts[gdatmodi.indxpoplmodi] or \
                                                     gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] != gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]):
        
        if rand() < gdat.probbrde or brth or deth:
            
            ## births and deaths
            if gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi] or deth:
                gdatmodi.propdeth = True
            elif gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.minmnumbpnts[gdatmodi.indxpoplmodi] or brth:
                gdatmodi.propbrth = True
            else:
                if rand() < 0.5:
                    gdatmodi.propbrth = True
                else:
                    gdatmodi.propdeth = True
            
            if gdatmodi.propbrth:
                gdatmodi.thisindxproptype = gdat.indxproptypebrth
            else:
                gdatmodi.thisindxproptype = gdat.indxproptypedeth

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
                gdatmodi.thisindxproptype = gdat.indxproptypesplt
            else:
                gdatmodi.thisindxproptype = gdat.indxproptypemerg
    else:
        if gdat.propwithsing:
            # temp
            gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal, gdatmodi.thisindxsampflux, gdatmodi.thisindxsampsind, \
                        gdatmodi.thisindxsampcurv, gdatmodi.thisindxsampexpo, gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxpntsfull, gdat.spectype) 
            gdatmodi.thisindxsampfull = concatenate((gdat.indxfixp, concatenate(gdatmodi.thisindxsampcomp)))
            gdatmodi.indxsampmodi = choice(gdatmodi.thisindxsampfull)

        gdatmodi.thisindxproptype = gdat.indxproptypewith
        gdatmodi.propwith = True

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

    oaxitype = len(fwhm.shape) == 3
    cntsbackfwhm = zeros_like(fwhm)
    for c in gdat.indxback:
        indxbacp = c * gdat.numbener + gdat.indxener
        if oaxitype:
            cntsback = bacp[indxbacp, None, None, None] * gdat.backflux[c][:, :, :, None] * gdat.expo[:, :, :, None] * \
                                                                                                gdat.deltener[:, None, None, None] * pi * fwhm[:, None, :, :]**2 / 4.
        else:
            cntsback = bacp[indxbacp, None, None] * gdat.backflux[c] * gdat.expo * pi * fwhm[:, None, :]**2 / 4.
            if gdat.enerdiff:
                cntsback *= gdat.deltener[:, None, None]
        cntsbackfwhm += mean(cntsback, 1)

    return cntsbackfwhm


def retr_sigm(gdat, cnts, cntsbackfwhm, lgal=None, bgal=None):
   
    oaxitype = len(cntsbackfwhm.shape) == 3
    if cnts.ndim == 2:
        if oaxitype:
            sigm = cnts / sum(cntsbackfwhm[:, :, 0], 1)[:, None]
        else:
            sigm = cnts / sum(cntsbackfwhm, 1)[:, None]
    else:
        if oaxitype:
            indxoaxitemp = retr_indxoaxipnts(gdat, lgal, bgal)
            sigm = cnts / swapaxes(cntsbackfwhm[:, :, indxoaxitemp], 1, 2)
        else:
            sigm = cnts / cntsbackfwhm[:, None, :]

    return sigm


def retr_probpois(data, modl):
    
    prob = data * log(modl) - modl - sp.special.gammaln(data + 1)

    return prob
    
        
def retr_sampvarb(gdat, indxpntsfull, samp, strg):
    
    if strg == 'true':
        strgtype = strg
    else:
        strgtype = ''
    spectype = getattr(gdat, strgtype + 'spectype')
    minmflux = getattr(gdat, strgtype + 'minmflux')
    binsflux = getattr(gdat, strgtype + 'binsflux')
    indxsamplgal, indxsampbgal, indxsampflux, indxsampsind, indxsampcurv, indxsampexpo, indxsampcomp = retr_indxsampcomp(gdat, indxpntsfull, spectype) 
    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxfixpnumbpnts] = samp[gdat.indxfixpnumbpnts]
    
    for k in gdat.indxfixp:
        sampvarb[k] = icdf_fixp(gdat, '', samp[k], k)
    
    if gdat.numbtrap > 0:
        for l in gdat.indxpopl:
            sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgang, 2. * gdat.maxmgang)
            sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgang, 2. * gdat.maxmgang) 
            
            fluxunit = samp[indxsampflux[l]]
            if gdat.fluxdisttype[l] == 'powr':
                sampvarb[indxsampflux[l]] = icdf_flux_powr(fluxunit, minmflux, gdat.maxmflux, sampvarb[gdat.indxfixpfluxdistslop[l]])
            if gdat.fluxdisttype[l] == 'bind':
                fluxdistnorm = sampvarb[gdat.indxfixpfluxdistnorm[l, :]]
                sampvarb[indxsampflux[l]] = icdf_bind(fluxunit, minmflux, gdat.maxmflux, binsflux, fluxdistnorm)
           
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
    gdat.exproaxitype = False


def retr_sdsspsfn(gdat):
   
    gdat.exprpsfp = array([0.25 / gdat.anglfact, 1.7e6, 1.9, 0.25 / gdat.anglfact, 2.1e6, 2.])
    gdat.exproaxitype = False


def retr_chanpsfn(gdat):

    gdat.exprpsfp = array([0.3 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.6e-1, 2.])
    gdat.exproaxitype = True
   

def retr_fermpsfn(gdat):
   
    gdat.exproaxitype = False
    
    gdat.fermreco = 7
    
    if gdat.fermreco == 8:
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
        if gdat.fermreco == 8:
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
    gdat.exprpsfp = zeros(gdat.numbener * gdat.numbevtt * numbpsfpform)
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
    
    # update the log-prior -- needed for diagnostics
    gdatmodi.thislpri = copy(gdatmodi.nextlpri)
        
    # update the log-prior sum
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
    #gdat.exprcnts = zeros((gdat.numbener, gdat.exprlgal.size, gdat.numbevtt))

    gdat.exprspec[0, 0, :] = fluxchansoft * 0.624e9
    gdat.exprspec[0, 1, :] = fluxchanhard * 0.624e9 / 16.
    
    # temp
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :]

    # temp
    gdat.exprspec[where(gdat.exprspec < 0.)] = 0.

    #gdat.exprcnts[0, :, 0] = cntschansoft
    #gdat.exprcnts[1, :, 0] = cntschanhard
        
    print 'gdat.exprspec[0, 0, :]'
    print gdat.exprspec[0, 0, :]
    print 'gdat.exprspec[0, 1, :]'
    print gdat.exprspec[0, 1, :]
    print 'log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :])'
    print log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :])
    print

    gdat.exprsind = -log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :]) / log(gdat.meanener[1] / gdat.meanener[0])
    gdat.exprsind[where(logical_not(isfinite(gdat.exprsind)))[0]] = 2.

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
    

    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.exprspec[where(isfinite(gdat.exprspec) == False)] = 0.
   
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3strg = fgl3['Source_Name']
    fgl3strgclss = fgl3['CLASS1']
    fgl3strgassc = fgl3['ASSOC1']
    
    fgl3timevari = fgl3['Variability_Index']

    fgl3spectype = fgl3['SpectrumType']
    gdat.exprsind = fgl3['Spectral_Index']
    gdat.exprcurv = fgl3['beta']
    gdat.exprexpo = fgl3['Cutoff'] * 1e-3
   
    #indxtimevari = where((fgl3timevari < 100.) & (gdat.exprspec[0, gdat.indxenerfluxdist[0], :] > gdat.minmflux))[0]
    indxtimevari = where(gdat.exprspec[0, gdat.indxenerfluxdist[0], :] > gdat.minmflux)[0]

    gdat.exprlgal = gdat.exprlgal[indxtimevari]
    gdat.exprbgal = gdat.exprbgal[indxtimevari]
    gdat.exprsind = gdat.exprsind[indxtimevari]
    gdat.exprcurv = gdat.exprcurv[indxtimevari]
    gdat.exprexpo = gdat.exprexpo[indxtimevari]
    gdat.exprspec = gdat.exprspec[:, :, indxtimevari]


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
        gdatmodi.nextsamp[indxsamp] = gdatmodi.nextsamp[indxsamp] % 1.
            
       
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

        dictmodi = {}
        dictmodi['indxsamplgal'], dictmodi['indxsampbgal'], dictmodi['indxsampflux'], dictmodi['indxsampsind'], \
        dictmodi['indxsampcurv'], dictmodi['indxsampexpo'], dictmodi['indxsampcomp'] = retr_indxsampcomp(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
        
        gdatmodi.indxparaprop = concatenate((gdat.indxfixpprop, concatenate(dictmodi['indxsampcomp'])))
        if k in concatenate(dictmodi['indxsamplgal']):
            print
        try:
            strgboolmodi = '%s' % (k in gdatmodi.indxsampmodi)
        except:
            strgboolmodi = ''
        print '%22s %14.4g %14.4g %14.4g %14.4g %14.4g %14s %10d' % (name, getattr(gdatmodi, 'thissamp')[k], getattr(gdatmodi, 'nextsamp')[k], gdatmodi.thissampvarb[k], \
                                               gdatmodi.nextsampvarb[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k], strgboolmodi, gdat.indxstdppara[k])
    

def retr_prop(gdat, gdatmodi, thisindxpnts=None):
 
    if gdat.verbtype > 1:
        print 'retr_prop()'
    
    gdatmodi.nextsamp = copy(gdatmodi.thissamp)
    gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
    gdatmodi.nextindxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
  
    if gdatmodi.propwith:
        
        gdatmodi.dictmodi = {}
        gdatmodi.dictmodi['indxsamplgal'], gdatmodi.dictmodi['indxsampbgal'], gdatmodi.dictmodi['indxsampflux'], gdatmodi.dictmodi['indxsampsind'], \
          gdatmodi.dictmodi['indxsampcurv'], gdatmodi.dictmodi['indxsampexpo'], gdatmodi.dictmodi['indxsampcomp'] = retr_indxsampcomp(gdat, gdatmodi.thisindxpntsfull, gdat.spectype)
        
        for k in gdat.indxprop:
            if gdatmodi.optipropdone:
                stdvstdp = gdatmodi.stdvstdp[k]
            else:
                stdvstdp = gdatmodi.nextstdvstdp[k]
            retr_gaus(gdat, gdatmodi, gdat.indxfixpprop[k], stdvstdp)
        
        for k in gdat.indxfixpdist:
            gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, '', gdatmodi.nextsamp[k], k)

        # rescale the unit sample vector due to the hyperparameter change
        if gdat.propwithsing:
            if gdatmodi.indxsampmodi in gdat.indxfixpdist:
                indxpoplsampmodi = [gdat.indxpoplsamp[gdatmodi.indxsampmodi]]
            else:
                indxpoplsampmodi = []
        else:
            indxpoplsampmodi = gdat.indxpopl
        
        for l in indxpoplsampmodi:
            
            ## flux distribution
            flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampflux[l]]
            if gdat.fluxdisttype[l] == 'powr':
                fluxunit = cdfn_flux_powr(flux, gdat.minmflux, gdat.maxmflux, gdatmodi.nextsampvarb[gdat.indxfixpfluxdistslop[l]])
            if gdat.fluxdisttype[l] == 'bind':
                fluxunit = cdfn_bind(flux, gdat.minmflux, gdat.maxmflux, gdat.binsflux, gdatmodi.nextsampvarb[gdat.indxfixpfluxdistnorm[l, :]])
            
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

        # PSs
        gdatmodi.thislfctprop = 0.
        for l in gdat.indxpopl:
            for k in range(len(gdatmodi.thisindxpntsfull[l])):
                for strg in gdat.liststrgcomp[l]:
                    if strg == 'flux':
                        gdatmodi.thisstdv = gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' +strg)] / \
                                                                                (gdatmodi.thissampvarb[gdatmodi.dictmodi['indxsampflux'][l][k]] / gdat.minmflux)**2.
                        gdatmodi.nextstdv = gdatmodi.stdvstdp[gdat.indxstdpflux] / (gdatmodi.nextsampvarb[gdatmodi.dictmodi['indxsampflux'][l][k]] / gdat.minmflux)**0.5
                        gdatmodi.thislfctprop += sum(0.5 * (gdatmodi.nextsamp[gdatmodi.dictmodi['indxsampflux'][l]] - \
                                               gdatmodi.thissamp[gdatmodi.dictmodi['indxsampflux'][l]])**2 * (1. / gdatmodi.thisstdv**2 - 1. / gdatmodi.nextstdv**2))
                    else:
                        gdatmodi.thisstdv = gdatmodi.stdvstdp[getattr(gdat, 'indxstdp' +strg)] / \
                                                                                (gdatmodi.thissampvarb[gdatmodi.dictmodi['indxsampflux'][l][k]] / gdat.minmflux)**0.5
                    retr_gaus(gdat, gdatmodi, gdatmodi.dictmodi['indxsamp' + strg][l][k], gdatmodi.thisstdv)
                        
        # temp
        gdatmodi.thislfctprop = 0.

    else:
        # number of unit sample vector elements to be modified
        numbcompmodi = gdat.numbcomp[gdatmodi.indxpoplmodi]
       
    if gdatmodi.propbrth:
       
        # find an empty slot in the PS list
        for k in range(gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]):
            if not k in gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi]:
                indxpntsbrth = k
                break
   
        # sample indices to add the new PS
        gdatmodi.indxsamptran = retr_indxsamppnts(gdat, gdatmodi.indxpoplmodi, array([indxpntsbrth]))
        
        # sample auxiliary variables
        gdatmodi.auxipara = rand(numbcompmodi)
        gdatmodi.nextsamp[gdatmodi.indxsamptran] = gdatmodi.auxipara
                
        # change the number of PS
        gdatmodi.nextsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].append(indxpntsbrth)
    
    # death
    if gdatmodi.propdeth:
        
        if thisindxpnts != None:
            dethindxindxpnts = thisindxpnts
        else:
            # occupied PS index to be killed
            dethindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
            
        # PS index to be killed
        indxpntsdeth = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][dethindxindxpnts]
    
        # sample indices to add the new PS
        gdatmodi.indxsamptran = gdat.indxsampcompinit + gdat.numbtrapcuml[gdatmodi.indxpoplmodi] + indxpntsdeth * gdat.numbcomp[gdatmodi.indxpoplmodi] + \
                                                                                                                                 gdat.indxcomp[gdatmodi.indxpoplmodi]
    
        # remove the PS from the occupied PS list
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].remove(indxpntsdeth)
    
        # change the number of PS
        gdatmodi.nextsamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
                  
    if gdat.numbtrap > 0:
        if gdatmodi.propwith:
            gdatmodi.indxsampmodi = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.dictmodi['indxsampcomp'])))
            gdatmodi.indxsampchec = gdatmodi.indxsampmodi
        else:
            gdatmodi.indxsampchec = []
            if gdatmodi.propbrth:
                gdatmodi.indxsampmodi = concatenate((gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran))
            else:
                gdatmodi.indxsampmodi = gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi, None]
    else:
        gdatmodi.indxsampchec = gdat.indxfixpprop
        gdatmodi.indxsampmodi = gdat.indxfixpprop
    
    # reject the sample if proposal is outside the prior
    indxchecfail = where((gdatmodi.nextsamp[gdatmodi.indxsampchec] < 0.) | (gdatmodi.nextsamp[gdatmodi.indxsampchec] > 1.))[0]
    if indxchecfail.size > 0:
        if gdat.verbtype > 1:
            for k in range(20):
                print 'Proposal rejected due to proposal outside the prior'
            print 'indxchecfail'
            print indxchecfail
            print
        gdatmodi.thisaccpprio = False
        print 'gdatmodi.thisaccpprio'

    if gdatmodi.thisaccpprio and not gdatmodi.propdeth:
        gdatmodi.nextsampvarb = retr_sampvarb(gdat, gdatmodi.nextindxpntsfull, gdatmodi.nextsamp, 'next')
   
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbpntsmodi = 3
        
        gdatmodi.nextsampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        gdatmodi.spltindxindxpnts = choice(thisindxindxpnts)
    
        # update the sample vector
        gdatmodi.indxsampfrst = gdat.indxsampcompinit + gdat.numbtrap * gdatmodi.indxpoplmodi + \
                                                int(gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]) * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp[gdatmodi.indxpoplmodi]

        gdatmodi.indxsampseco = gdat.indxsampcompinit + gdat.numbtrap * gdatmodi.indxpoplmodi + indxpntsbrth * gdat.numbcomp[gdatmodi.indxpoplmodi]
        indxfinlseco = gdatmodi.indxsampseco + gdat.numbcomp[gdatmodi.indxpoplmodi]
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, indxfinlfrst, dtype=int), arange(gdatmodi.indxsampseco, indxfinlseco, dtype=int)))
        
        thislgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisbgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, gdatmodi.spltindxindxpnts]]
        gdatmodi.fluxpare = thisspec[gdat.indxenerfluxdist[0]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, gdatmodi.spltindxindxpnts]]
        gdatmodi.fluxpare = thisspec[gdat.indxenerfluxdist[0]]
        thisspep = gdatmodi.thissampvarb[gdatmodi.thisindxsampspep[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts, :]]
        
        # determine the new element parameters
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
            gdatmodi.thisaccpprio = False

        # calculate the list of pairs
        ## current
        gdatmodi.thislistpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], \
                                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]])
        gdatmodi.thisnumbpair = len(gdatmodi.thislistpair)
        
        if gdatmodi.thisaccpprio:

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

            ## first new element
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcomplgal] = cdfn_self(gdatmodi.spltlgalfrst, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompbgal] = cdfn_self(gdatmodi.spltbgalfrst, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompflux] = cdfn_flux_powr(gdatmodi.fluxfrst, gdat.minmflux, gdat.maxmflux, \
                                                                                gdatmodi.thissampvarb[gdat.indxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompsind] = cdfn_gaus(gdatmodi.spltsindfrst, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
    
            nextspecfrst = retr_spec(gdat, gdatmodi.fluxfrst, spep=gdatmodi.spltsindfrst, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

            ## second new element
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
            gdatmodi.indxsampfrst = gdat.indxsampcompinit + gdat.numbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxfrst
            indxfinlfrst = gdatmodi.indxsampfrst + gdat.numbcomp[gdatmodi.indxpoplmodi]
            
            ## second PS
            gdatmodi.indxsampseco = gdat.indxsampcompinit + gdat.numbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxseco
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
        
    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        
        ## Jacobian
        jcbnfacttemp = log(gdatmodi.fluxpare * fabs(gdatmodi.auxipara[1] * (sin(gdatmodi.auxipara[2]) * cos(gdatmodi.auxipara[2]) + cos(gdatmodi.auxipara[2])**2)))
        if gdatmodi.propsplt:
            gdatmodi.thisjcbnfact = jcbnfacttemp
        else:
            gdatmodi.thisjcbnfact = -jcbnfacttemp
        
        ## combinatorial factor
        thisnumbpnts = gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]
        if gdatmodi.propsplt:
            gdatmodi.thiscombfact = log(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2 / gdatmodi.nextnumbpair)
        else:
            gdatmodi.thiscombfact = log(gdatmodi.thisnumbpair / gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2)

        
def retr_indxsamppnts(gdat, l, indxpnts):

    indxsamppnts = gdat.indxsampcompinit + gdat.numbtrapcuml[l] + indxpnts[None, :] * gdat.numbcomp[l] + gdat.indxcomp[l][:, None]
    indxsamppnts = indxsamppnts.flatten()

    return indxsamppnts


def retr_factoaxi(gdat, bins, norm, indx):

    factoaxi = 1. + norm[:, None, None] * (bins[None, None, :] / gdat.oaxipivt)**indx[:, None, None]
    
    return factoaxi


def gang_detr():

    gang, aang, lgal, bgal = sympy.symbols('gang aang lgal bgal')

    AB = sympy.matrices.Matrix([[a1*b1,a1*b2,a1*b3],[a2*b1,a2*b2,a2*b3],[a3*b1,a3*b2,a3*b3]])
    print AB.det()


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, binsoaxi=None, oaxitype=None, strgpara=''):

    numbpsfpform = getattr(gdat, strgpara + 'numbpsfpform')
    numbpsfptotl = getattr(gdat, strgpara + 'numbpsfptotl')
    
    indxpsfpinit = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    if oaxitype:
        indxpsfponor = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp]
        indxpsfpoind = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp] + 1
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * arcsin(sqrt(2. - 2. * cos(gdat.binsanglplot)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        if oaxitype:
            scalangl = thisangl[None, :, None, None]
        else:
            scalangl = thisangl[None, :, None]
    
    if oaxitype:
        factoaxi = retr_factoaxi(gdat, binsoaxi, psfp[indxpsfponor], psfp[indxpsfpoind])
   
    if psfntype == 'singgaus':
        sigc = psfp[indxpsfpinit]
        if oaxitype:
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
        else:
            sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
    
    elif psfntype == 'singking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        psfn = retr_singking(scalangl, sigc, gamc)
        if oaxitype:
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
            gamc = gamc[:, None, :, None]
        else:
            sigc = sigc[:, None, :]
            gamc = gamc[:, None, :]
        
    elif psfntype == 'doubgaus':
        frac = psfp[indxpsfpinit]
        sigc = psfp[indxpsfpinit+1]
        sigt = psfp[indxpsfpinit+2]
        if oaxitype:
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
        if oaxitype:
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
        if oaxitype:
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
        fact = 2. * pi * trapz(psfnnorm * sin(gdat.binsanglplot[None, :, None]), gdat.binsanglplot, axis=1)[:, None, :]
        psfn /= fact

    # temp
    if True and (gdat.strgcnfg == 'pcat_ferm_expr_ngal' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp1' or \
                                                            gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp2' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp3'):
        print 'CORRECTING THE PSF.'
        tempcorr = array([1., 0.8, 0.8])
        psfn *= tempcorr[:, None, None]

    return psfn


def retr_axis(gdat, strgvarb, minm=None, maxm=None, numb=None, bins=None, scal='self', invr=False, strg=''):
    
    if bins == None:
        if scal == 'self':
            bins = linspace(minm, maxm, numb + 1)
            meanvarb = (bins[1:] + bins[:-1]) / 2.
        else:
            bins = logspace(log10(minm), log10(maxm), numb + 1)
            meanvarb = sqrt(bins[1:] * bins[:-1])
    else:
        numb = bins.size - 1

    if invr:
        bins = bins[::-1]

    if scal == 'self':
        meanvarb = (bins[1:] + bins[:-1]) / 2.
    else:
        meanvarb = sqrt(bins[1:] * bins[:-1])
    
    indx = arange(numb)
    delt = diff(bins) 

    setattr(gdat, strg + 'bins' + strgvarb, bins)
    setattr(gdat, strg + 'mean' + strgvarb, meanvarb)
    setattr(gdat, strg + 'delt' + strgvarb, delt)
    setattr(gdat, strg + 'numb' + strgvarb, numb)
    setattr(gdat, strg + 'indx' + strgvarb, indx)


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

    psecodim = zeros(gdat.numbsidecart / 2)
    for k in gdat.indxwvecodimplot:
        indxwvec = where((gdat.meanwvec > gdat.binswvecodimplot[k]) & (gdat.meanwvec < gdat.binswvecodimplot[k+1]))
        psecodim[k] = mean(psec[indxwvec])
    
    psecodim *= gdat.meanwvecodimplot**2

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


def retr_eerrnorm(minmvarb, maxmvarb, meanvarb, stdvvarb):
   
    cdfnminm = 0.5 * (sp.special.erf((minmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfnmaxm = 0.5 * (sp.special.erf((maxmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfndiff = cdfnmaxm - cdfnminm
    
    return cdfnminm, cdfndiff
    

def retr_detrcatl(gdat):
  
    if gdat.verbtype > 1:
        print 'gdat.listindxpntsfull'
        print gdat.listindxpntsfull

    # find the total number of elements across all the samples
    numbelem = 0
    indxtupl = []
    indxelem = []
    indxelemsamp = []
    for n in gdat.indxsamptotl:
        indxelem.append([])
        indxelemsamptemp = []
        for l in gdat.indxpopl:
            indxelem[n].append([])
            for k in range(len(gdat.listindxpntsfull[n][l])):
                indxelem[n][l].append(numbelem)
                indxelemsamptemp.append(numbelem)
                indxtupl.append([n, l, k])
                numbelem += 1
        indxelemsamp.append(array(indxelemsamptemp))
    
    if gdat.verbtype > 1:
        print 'indxelem'
        print indxelem
        print 'indxtupl'
        print indxtupl
        print 'indxelemsamp'
        print indxelemsamp
        print 'numbelem'
        print numbelem

    cntr = -1
    numbassc = zeros(numbelem)
    listdist = zeros((numbelem, numbelem)) + 1e20
    cntrpercsave = 0
    for n in gdat.indxsamptotl:
        for l in gdat.indxpopl:
            for k in range(len(gdat.listindxpntsfull[n][l])):
                indxelemfrst = indxelem[n][l][k]
                for nn in gdat.indxsamptotl:
                    for ll in gdat.indxpopl:
                        indxpntsseco = array(gdat.listindxpntsfull[nn][ll])
                        indxelemseco = array(indxelem[nn][ll])
                        
                        cntr += indxpntsseco.size
                        
                        if n == nn or len(gdat.listindxpntsfull[nn][ll]) == 0:
                            continue
                
                        ## compute distance
                        indxsampfrst = retr_indxsamppnts(gdat, l, array([gdat.listindxpntsfull[n][l][k]]))
                        indxsampseco = retr_indxsamppnts(gdat, l, indxpntsseco)
                        dist = sqrt(sum((gdat.listsamp[n, indxsampfrst][:, None] - gdat.listsamp[nn, indxsampseco.reshape((indxsampfrst.size, -1))])**2, axis=0))
                        
                        listdist[indxelemfrst, indxelemseco] = dist
                        
                        if gdat.verbtype > 1:
                            if cntr % 10000 == 0:
                                print 'cntr'
                                print cntr

                        cntrperc = 5 * floor(20. * cntr / numbelem**2)
                        if cntrperc > cntrpercsave:
                            if gdat.verbtype > 0:
                                print 'Distance table calculation %d%% completed.' % cntrpercsave
                            cntrpercsave = cntrperc

    indxelemleft = range(numbelem)
    distthrs = 0.1

    # list of sample lists of the labeled element
    indxelemassc = []
    cntr = 0
    
    if gdat.verbtype > 1:
        print 'listdist'
        print listdist
    
    while len(indxelemleft) > 0 and cntr < 10:
        
        if gdat.verbtype > 1:
            print 'indxelemassc'
            print indxelemassc

        # count number of associations
        numbassc = zeros(numbelem, dtype=int) - 1
        for p in range(len(indxelemleft)):
            numbassc[indxelemleft[p]] = where(listdist[indxelemleft[p], indxelemleft] < distthrs)[0].size
        
        # determine the element with the highest number of neighbors
        indxelemcntr = argmax(numbassc)
        indxsamptotlcntr = indxtupl[indxelemcntr][0]
        indxpoplcntr = indxtupl[indxelemcntr][1]
        indxpntscntr = indxtupl[indxelemcntr][2]

        indxelemassc.append([])
        indxelemassc[cntr].append(indxelemcntr)
        indxelemleft.remove(indxelemcntr)

        if gdat.verbtype > 1:
            print 'Match step %d' % cntr
            print 'numbassc'
            print numbassc
            print 'indxelemcntr'
            print indxelemcntr
            print 'indxelemleft'
            print indxelemleft
        
        # add the central element sample
        # add the associated element samples
        if len(indxelemleft) > 0:
            for n in gdat.indxsamptotl:
                
                indxelemtemp = intersect1d(array(indxelemleft), indxelemsamp[n])
                
                if gdat.verbtype > 1:
                    print 'n'
                    print n
                    print 'indxelemsamp[n]'
                    print indxelemsamp[n]
                    print 'indxelemtemp'
                    print indxelemtemp
                
                if n == indxsamptotlcntr:
                    continue
                
                if indxelemtemp.size > 0:
                    indxleft = argsort(listdist[indxelemcntr, indxelemtemp])[0]
                    indxelemthis = indxelemtemp[indxleft]
                
                    if gdat.verbtype > 1:
                        print 'indxleft'
                        print indxleft
                        print 'indxelemthis'
                        print indxelemthis
                        print 'listdist[indxelemcntr, indxelemthis]'
                        print listdist[indxelemcntr, indxelemthis]
                
                    if listdist[indxelemcntr, indxelemthis] < distthrs:
                        indxelemassc[cntr].append(indxelemthis)
                        indxelemleft.remove(indxelemthis)
                        if gdat.verbtype > 1:
                            print 'Appending...'
            
                if gdat.verbtype > 1:
                    print 

            cntr += 1
        
        if gdat.verbtype > 1:
            print 
            print 
            print 
        
    if gdat.verbtype > 1:
        print 'gdat.listindxpntsfull'
        print gdat.listindxpntsfull
        print 'indxelemassc'
        print indxelemassc
        print
        for strgcomp in gdat.liststrgcomptotl:
            print 'strgcomp'
            print strgcomp
            print 'getattr(gdat, list + strgcomp)'
            print getattr(gdat, 'list' + strgcomp)

    gdat.dictglob['listelemdetr'] = []
    gdat.dictglob['postelemdetr'] = []
    for r in range(len(indxelemassc)): 
        gdat.dictglob['listelemdetr'].append([])
        gdat.dictglob['listelemdetr'][r] = {}
        gdat.dictglob['postelemdetr'].append([])
        gdat.dictglob['postelemdetr'][r] = {}
        for strgcomp in gdat.liststrgcomptotl:
            gdat.dictglob['listelemdetr'][r][strgcomp] = []
            gdat.dictglob['postelemdetr'][r][strgcomp] = zeros(3)
        for k in range(len(indxelemassc[r])):
            indxsamptotlcntr = indxtupl[indxelemassc[r][k]][0]
            indxpoplcntr = indxtupl[indxelemassc[r][k]][1]
            indxpntscntr = indxtupl[indxelemassc[r][k]][2]
            
            if gdat.verbtype > 1:
                print 'indxsamptotlcntr'
                print indxsamptotlcntr
                print 'indxpoplcntr'
                print indxpoplcntr
                print 'indxpntscntr'
                print indxpntscntr
                print
            
            for strgcomp in gdat.liststrgcomptotl:
                temp = getattr(gdat, 'list' + strgcomp)
                temp = temp[indxpoplcntr][indxsamptotlcntr][indxpntscntr]
                gdat.dictglob['listelemdetr'][r][strgcomp].append(temp)
    
    gdat.numbpntsdetr = len(gdat.dictglob['listelemdetr'])

    for r in range(len(indxelemassc)): 
        for strgcomp in gdat.liststrgcomptotl:
            arry = array(gdat.dictglob['listelemdetr'][r][strgcomp])
            gdat.dictglob['postelemdetr'][r][strgcomp][0] = median(arry)
            gdat.dictglob['postelemdetr'][r][strgcomp][1] = percentile(arry, 16.)
            gdat.dictglob['postelemdetr'][r][strgcomp][2] = percentile(arry, 84.)
   
    if gdat.verbtype > 1:
        print 'gdat.dictglob[listelemdetr]'
        for r in range(len(gdat.dictglob['listelemdetr'])):
            print 'r'
            print r
            for strgcomp in gdat.liststrgcomptotl:
                print strgcomp
                print gdat.dictglob['listelemdetr'][r][strgcomp]
                print 
        print 'gdat.dictglob[postelemdetr]'
        print gdat.dictglob['postelemdetr']
        
        
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


def setp_indxswepsave(gdat):

    gdat.indxswep = arange(gdat.numbswep)
    gdat.boolsave = zeros(gdat.numbswep, dtype=bool)
    gdat.indxswepsave = arange(gdat.numbburn, gdat.numbburn + gdat.numbsamp * gdat.factthin, gdat.factthin)
    gdat.boolsave[gdat.indxswepsave] = True
    gdat.indxsampsave = zeros(gdat.numbswep, dtype=int) - 1
    gdat.indxsampsave[gdat.indxswepsave] = arange(gdat.numbsamp)
    

def retr_pntscnts(gdat, lgal, bgal, spec):
    
    indxpixltemp = retr_indxpixl(gdat, bgal, lgal)
    cnts = zeros((gdat.numbener, lgal.size, gdat.numbevtt))
    for k in range(lgal.size):
        cnts[:, k, :] += spec[:, k, None] * gdat.expo[:, indxpixltemp[k], :]
    if gdat.enerdiff:
        cnts *= gdat.deltener[:, None, None]
    
    return cnts


def retr_critmden(gdat, adissour, adishost, adishostsour):
    
    critmden = gdat.factnewtlght / 4. / pi * adissour / adishostsour / adishost
        
    return critmden


def retr_massfrombein(gdat, adissour, adishost, adishostsour):

    critmden = retr_critmden(gdat, adissour, adishost, adishostsour)
    massfrombein = pi * adishost**2 * critmden

    return massfrombein


def retr_massfromdeflscal(gdat, adissour, adishost, adishostsour, anglscal, anglcutf):
    
    critmden = retr_critmden(gdat, adissour, adishost, adishostsour)
    
    fraccutf = anglcutf / anglscal
    massfromdeflscal = pi * anglscal * adishost**2 * critmden / (1. - log(2.)) * fraccutf**2 / (fraccutf**2 + 1.)**2 * \
                                                                                ((fraccutf**2 - 1.) * log(fraccutf) + fraccutf * pi - (fraccutf**2 + 1.))
        
    return massfromdeflscal


def setpinit(gdat, boolinitsetp=False):

    if gdat.pntstype == 'lens' and gdat.strgproc == 'fink2.rc.fas.harvard.edu':
        cliblens = ctypes.CDLL(os.environ["PCAT_PATH"] + '/cliblens.so')
        cliblens.retr_deflsubh()

    # samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    
    # samples to be saved from all chains
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc
    
    # run tag
    gdat.rtag = retr_rtag(gdat)
    
    gdat.liststrgtype = ['spatdisttype', 'fluxdisttype', 'spectype']

    # paths
    ## data
    gdat.pathpixlcnvt = gdat.pathdata + 'pixlcnvt/'
    gdat.pathprox = gdat.pathdata + 'prox/'
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
        for strg in ['llik']:
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
    gdat.minmdatacnts = 0.
    
    ## number of bins in histogram plots
    gdat.numbbinsplot = 20
    ## number of bins in hyperprior plots
    gdat.numbbinsplotprio = 100
    # temp
    ### number of standard deviations away from the mean of spectral parameters to plot
    gdat.numbstdvspepdist = 3.
    ### minima and maxima for spectral parameters
    gdat.numbstdv = 3.
    
    retr_axis(gdat, 'flux', gdat.minmflux, gdat.maxmflux, gdat.numbfluxdistnorm - 1, scal='logt')

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
    gdat.labldeltllik = r'\Delta_a P(D|x)'
    gdat.labldeltllikunit = ''
    gdat.labldistsour = 'd_{sa}'
    gdat.labldistsourunit = gdat.lablgangunit
    gdat.labldotpsour = r'\vec{\alpha}_a \cdot \nabla f_s'
    gdat.labldotpsourunit = ''
    
    ## descriptive string for the elements 
    if gdat.pntstype == 'lght':
        gdat.strgelem = 'PS'
    if gdat.pntstype == 'lens':
        gdat.strgelem = 'Subhalo'
    
    # scalar variables
    gdat.minmfracsubh = 0.
    gdat.maxmfracsubh = 0.3
    gdat.scalfracsubh = 'self'

    # set up the indices of the fitting model
    retr_indxsamp(gdat)

    # construct the fitting model
    setp_fixp(gdat)
    
    # for each parameter in the fitting model, determine if there is a corresponding parameter in the generative model
    gdat.corrfixp = empty(gdat.numbfixp)
    for k in gdat.indxfixp:
        try:
            gdat.corrfixp[k] = getattr(gdat, 'true' + gdat.namefixp[k])
        except:
            gdat.corrfixp[k] = None
   
    # copy fixp variables to individual variables
    for k, namefixp in enumerate(gdat.namefixp):
        setattr(gdat, 'minm' + namefixp, gdat.minmfixp[k])
        setattr(gdat, 'maxm' + namefixp, gdat.maxmfixp[k])
        setattr(gdat, 'scal' + namefixp, gdat.scalfixp[k])
    
    for l in gdat.indxpopl:
        setattr(gdat, 'minmnumbpntspop%d' % l, gdat.minmnumbpnts[l])
        setattr(gdat, 'maxmnumbpntspop%d' % l, gdat.maxmnumbpnts[l])
    
    # list of scalar variable names
    gdat.liststrgvarbscal = list(gdat.namefixp)
    if gdat.pntstype == 'lens':
        gdat.liststrgvarbscal += ['fracsubh']

    # construct bins for the scalar variables
    for strgvarbscal in gdat.liststrgvarbscal:
        minm = getattr(gdat, 'minm' + strgvarbscal)
        maxm = getattr(gdat, 'maxm' + strgvarbscal)
        scal = getattr(gdat, 'scal' + strgvarbscal)
        if scal == 'logt' or scal == 'lgau':
            thisscal = 'logt'
        else:
            thisscal = 'self'
        retr_axis(gdat, strgvarbscal, minm, maxm, gdat.numbbinsplot, scal=thisscal)

    # temp
    gdat.anglscal = 0.05 / gdat.anglfact 
    gdat.anglcutf = 1. / gdat.anglfact 
    
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
        if strgfeat == 'flux' and gdat.pntstype == 'lens' or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal' or strgfeat == 'distsour':
            gdat.dictglob['fact' + strgfeat + 'plot'] = gdat.anglfact
        else:
            gdat.dictglob['fact' + strgfeat + 'plot'] = 1.
        setattr(gdat, 'numb' + strgfeat + 'plot', 20)
        
        if strgfeat == 'flux' or strgfeat == 'expo' or strgfeat == 'cnts':
            gdat.dictglob['scal' + strgfeat + 'plot'] = 'logt'
        else:
            gdat.dictglob['scal' + strgfeat + 'plot'] = 'self'

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

    if gdat.pntstype == 'lens':
        h5f = h5py.File(gdat.pathdata + 'inpt/adis.h5','r')
        reds = h5f['reds'][()]
        adis = h5f['adis'][()]
        adistdim = h5f['adistdim'][()]
        gdat.adisobjt = interp1d_pick(reds, adis)
        h5f.close()

        gdat.redshost = 0.5
        gdat.redssour = 2.
        gdat.adishost = gdat.adisobjt(gdat.redshost) * 1e3 # [kpc]
        gdat.adissour = gdat.adisobjt(gdat.redssour) * 1e3 # [kpc]
        gdat.adishostsour = gdat.adissour - (1. + gdat.redshost) / (1. + gdat.redssour) * gdat.adishost
        gdat.adisfact = gdat.adishost * gdat.adissour / gdat.adishostsour
        gdat.factnewtlght = 2.09e16 # Msun / kpc
        gdat.massfrombein = retr_massfrombein(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        gdat.critmden = retr_critmden(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        gdat.massfromdeflscal = retr_massfromdeflscal(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour, gdat.anglscal, gdat.anglcutf)
        gdat.minmmass = gdat.massfromdeflscal * gdat.minmflux
        gdat.maxmmass = gdat.massfromdeflscal * gdat.maxmflux
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
    
    # angular deviation
    # check if 1000 is too much
    # temp
    gdat.numbangl = 100
    if gdat.exprtype == 'ferm':
        gdat.maxmangl = 17. / gdat.anglfact
    if gdat.exprtype == 'chan':
        gdat.maxmangl = 15. / gdat.anglfact
    if gdat.exprtype == 'hubb':
        gdat.maxmangl = 2. / gdat.anglfact
    retr_axis(gdat, 'anglplot', 0., gdat.maxmangl, gdat.numbangl)
    retr_axis(gdat, 'anglfull', 0., 3. * gdat.maxmgang, gdat.numbangl)
    gdat.binsanglcosi = sort(cos(gdat.binsanglplot))
    
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
    
    # PS indices in each population
    gdat.indxpntspopl = []
    for l in gdat.indxpopl:
        gdat.indxpntspopl.append(arange(sum(gdat.maxmnumbpnts[:l]), sum(gdat.maxmnumbpnts[:l+1])))
    
    ## PSF class indices for which images will be plotted
    if gdat.numbevtt == 1:
        gdat.indxevttplot = gdat.indxevtt
    else:
        gdat.indxevttplot = concatenate((array([-1]), gdat.indxevtt))
    
    gdat.numbenerevtt = gdat.numbener * gdat.numbevtt
    
    gdat.correxpo = True
    
    # off-axis angle
    gdat.numboaxi = 10
    gdat.minmoaxi = 0.
    gdat.maxmoaxi = 1.1 * sqrt(2.) * gdat.maxmgang
    retr_axis(gdat, 'oaxiplot', gdat.minmoaxi, gdat.maxmoaxi, gdat.numboaxi)
    gdat.binsoaxiopen = gdat.binsoaxiplot[:-1]
    gdat.indxoaxi = arange(gdat.numboaxi)

    # pixelization
    if gdat.datatype == 'mock':
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2

    gdat.numbchrototl = 8
    gdat.numbchroproc = 3

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
    ## number of deflection components to plot
    gdat.numbdeflpntsplot = 3
    gdat.numbdeflsingplot = gdat.numbdeflpntsplot + 2

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
    
    gdat.binslgalcart = linspace(gdat.minmlgaldata, gdat.maxmlgaldata, gdat.numbsidecart + 1)
    gdat.binsbgalcart = linspace(gdat.minmbgaldata, gdat.maxmbgaldata, gdat.numbsidecart + 1)
    gdat.binslgalcartmesh, gdat.binsbgalcartmesh = meshgrid(gdat.binslgalcart, gdat.binsbgalcart)
    gdat.lgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
    gdat.bgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
    if gdat.pixltype == 'cart':
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
    
    # angular diameter distance
    if gdat.pntstype == 'lens':
        gdat.factkpcs = 1e3
        gdat.adis = gdat.adishost * gdat.factkpcs
    else:
        gdat.adis = None
    
    # power spectra
    if gdat.pixltype == 'cart':
        gdat.numbwvecodim = gdat.numbsidecart
        gdat.minmanglodim = 0.
        gdat.maxmanglodim = 2. * gdat.maxmgang
        gdat.minmmpolodim = 0.#1. / gdat.numbsidecart / gdat.sizepixl
        gdat.maxmmpolodim = 1. / 2. / gdat.sizepixl
        retr_axis(gdat, 'anglodimplot', gdat.minmanglodim, gdat.maxmanglodim, gdat.numbsidecart, invr=True)
        retr_axis(gdat, 'mpolodimplot', gdat.minmmpolodim, gdat.maxmmpolodim, gdat.numbsidecart / 2)
        if gdat.adis != None:
            gdat.minmwvecodim = gdat.minmmpolodim / gdat.adis
            gdat.maxmwvecodim = gdat.maxmmpolodim / gdat.adis
            gdat.minmwlenodim = gdat.minmanglodim * gdat.adis
            gdat.maxmwlenodim = gdat.maxmanglodim * gdat.adis
            retr_axis(gdat, 'wvecodimplot', gdat.minmwvecodim, gdat.maxmwvecodim, gdat.numbsidecart / 2)
            retr_axis(gdat, 'wlenodimplot', gdat.minmwlenodim, gdat.maxmwlenodim, gdat.numbsidecart, invr=True)
            gdat.meanwveclgal, gdat.meanwvecbgal = meshgrid(gdat.meanwvecodimplot, gdat.meanwvecodimplot, indexing='ij')
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
            if gdat.datatype == 'inpt':
                backfluxtemp = zeros_like(gdat.exprdataflux) + gdat.back[c]
        else:
            path = gdat.pathinpt + gdat.back[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'cart':
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        gdat.backflux.append(backfluxtemp)
    
    # only include desired energy and PSF class bins 
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    if gdat.correxpo:
        gdat.expo = gdat.expo[gdat.indxcubeincl]
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcubeincl]
  
    if gdat.pixltype == 'cart':
        gdat.backfluxcart = []
        for c in gdat.indxback:
            gdat.backfluxcart.append(gdat.backflux[c].reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt)))
        gdat.expocart = gdat.expo.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))

    for c in gdat.indxback:
        if amin(gdat.backflux[c]) <= 0.:
            raise Exception('Background templates must be positive.')

    if gdat.datatype == 'inpt':
        gdat.exprdataflux = gdat.exprdataflux[gdat.indxcubeincl]

    # temp
    if False and gdat.datatype == 'inpt':

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

    # exclude voxels with vanishing exposure
    if gdat.correxpo:
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
    
    if False and gdat.evalcirc != 'full' and gdat.numbpixl * gdat.maxmnumbpntstotl < 1e5:
        gdat.calcerrr = True
    else:
        gdat.calcerrr = False
    
    if gdat.pntstype == 'lght':
        gdat.exprpsfn = retr_psfn(gdat, gdat.exprpsfp, gdat.indxener, gdat.binsanglplot, gdat.exprpsfntype, gdat.binsoaxiplot, gdat.exproaxitype)
    
    if gdat.evalcirc != 'full':
        
        if gdat.evalcirc == 'psfn':
            gdat.numbfluxprox = 3
        if gdat.evalcirc == 'bein':
            gdat.numbfluxprox = 3
        
        gdat.indxprox = arange(gdat.numbfluxprox)
        gdat.binsprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
        
        if gdat.evalcirc == 'psfn':
            # determine the maximum angle at which the PS flux map will be computed
            gdat.maxmangleval = empty(gdat.numbfluxprox)
            for h in gdat.indxprox:
                if gdat.specfraceval == 0:
                    gdat.maxmangleval[h] = 3. * gdat.maxmgang
                else:  
                    frac = min(1e-2, gdat.specfraceval * gdat.binsprox[0] / gdat.binsprox[h+1])
                    psfnwdth = retr_psfnwdth(gdat, gdat.exprpsfn, frac)
                    gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                    gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
        
        if gdat.evalcirc == 'bein':
            gdat.maxmangleval = 10. * gdat.binsprox[1:]
            #gdat.maxmangleval = array([1. / gdat.anglfact]) # 4 * gdat.binsprox[1:]

        if gdat.pntstype == 'lght' and gdat.maxmangl - amax(gdat.maxmangleval) < 1.1 * sqrt(2) * (gdat.maxmgang - gdat.maxmgangdata):
            print 'gdat.maxmangl'
            print gdat.maxmangl * gdat.anglfact
            print 'gdat.maxmangleval'
            print gdat.maxmangleval * gdat.anglfact
            print 'gdat.maxmgang'
            print gdat.maxmgang * gdat.anglfact
            print 'gdat.maxmgangdata'
            print gdat.maxmgangdata * gdat.anglfact
            raise Exception('Angular axis is too short.')

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathprox + 'indxprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), \
                                                                                                            1e2 * amax(gdat.maxmangleval), gdat.numbfluxprox)
    
    if gdat.evalcirc != 'full':

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathprox + 'indxprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), \
                                                                                                            1e2 * amax(gdat.maxmangleval), gdat.numbfluxprox)
        if gdat.verbtype > 0 and boolinitsetp:
            print 'Element evaluation will be performed up to'
            print gdat.maxmangleval * gdat.anglfact

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
                    if indxpixlproxtemp.size > 1e4:
                        indxpixlproxtemp = -1
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
    gdat.alphpnts = 0.7
    gdat.alphmaps = 1.
    
    # number of colorbar ticks in the maps
    gdat.numbtickcbar = 11
    
    # turn off relevant proposal types
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

    gdat.indxstdplgal = gdat.numbfixpprop
    gdat.indxstdpbgal = gdat.numbfixpprop + 1
    gdat.indxstdpflux = gdat.numbfixpprop + 2
    if gdat.numbener > 1:
        gdat.indxstdpsind = gdat.numbfixpprop + 3
        gdat.indxstdpcurv = gdat.numbfixpprop + 4
        gdat.indxstdpexpo = gdat.numbfixpprop + 4
    gdat.numbstdp = gdat.numbfixpprop + gdat.maxmnumbcomp
    gdat.strgstdp = concatenate((array(gdat.strgfixp)[gdat.indxfixpprop], gdat.liststrgcomptotl))
    gdat.strgstdp = list(gdat.strgstdp)
    gdat.namestdp = concatenate((array(gdat.namefixp)[gdat.indxfixpprop], gdat.liststrgcomptotl))
    gdat.indxstdp = arange(gdat.numbstdp)
    gdat.indxstdpfixp = arange(gdat.numbfixpprop)
    gdat.indxstdpcomp = setdiff1d(gdat.indxstdp, gdat.indxstdpfixp)
    
    # proposal scale
    if gdat.pntstype == 'lens':
        gdat.stdvstdp = 1e-4 + zeros(gdat.numbstdp)
        gdat.stdvstdp[gdat.indxfixphypr+gdat.numbpopl] = 1e-3
        gdat.stdvstdp[gdat.indxstdpcomp] = 1e-3
    else:
        if gdat.exprtype == 'ferm':
            gdat.stdvstdp = 1e-3 + zeros(gdat.numbstdp)
            gdat.stdvstdp[gdat.indxfixphypr+gdat.numbpopl] = 1e-3
            gdat.stdvstdp[gdat.indxstdpcomp] = 1e-3
        if gdat.exprtype == 'chan':
            gdat.stdvstdp = 1e-3 + zeros(gdat.numbstdp)
            gdat.stdvstdp[gdat.indxfixphypr+gdat.numbpopl] = 1e-3
            gdat.stdvstdp[gdat.indxstdpcomp] = 1e-3

    # proposal scale indices for each parameter
    indxpntsfull = [range(gdat.maxmnumbpnts[l]) for l in gdat.indxpopl]
    gdat.indxsamplgal, gdat.indxsampbgal, gdat.indxsampflux, gdat.indxsampsind, gdat.indxsampcurv, gdat.indxsampexpo, \
                                                                                            gdat.indxsampcomp = retr_indxsampcomp(gdat, indxpntsfull, gdat.spectype) 
    gdat.indxstdppara = zeros(gdat.numbpara, dtype=int) - 1
    cntr = 0
    gdat.indxparaprop = zeros(gdat.numbfixpprop, dtype=int)
    for k in gdat.indxpara:
        if k in gdat.indxfixpprop:
            gdat.indxstdppara[k] = cntr
            gdat.indxparaprop[cntr] = k
            cntr += 1
        for l in gdat.indxpopl:
            for strgcomp in gdat.liststrgcomp[l]:
                if k in getattr(gdat, 'indxsamp' + strgcomp)[l]:
                    gdat.indxstdppara[k] = getattr(gdat, 'indxstdp' + strgcomp)


    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = retr_cntsmaps(gdat, gdat.exprdataflux)
    

    if gdat.pntstype == 'lght':
        gdat.limsangl = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        gdat.limspsfn = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.trueoaxitype:
                    psfn = gdat.exprpsfn[i, :, m, 0]
                else:
                    psfn = gdat.exprpsfn[i, :, m]
                maxmpsfn = amax(psfn)
                gdat.limsangl[i][m] = [0., gdat.binsanglplot[amax(where(psfn > 1e-6 * maxmpsfn)[0])] * gdat.anglfact]
                gdat.limspsfn[i][m] = [maxmpsfn * 1e-6, maxmpsfn]
       
    # spatially averaged background flux 
    gdat.backfluxmean = zeros((gdat.numbback, gdat.numbener))
    for c in gdat.indxback:
        for i in gdat.indxener:
            if gdat.correxpo:
                gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, :] * gdat.expo[i, :, :]) / sum(gdat.expo[i, :, :])
            else:
                gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, 0])
    
    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(10000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)
    gdat.indxcubesave = meshgrid(gdat.indxener, gdat.indxpixlsave, gdat.indxevtt, indexing='ij')
   
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
                print '%20s%25s%5s%20.6g%20.6g' % (gdat.namefixp[k], gdat.strgfixp[k], gdat.scalfixp[k], gdat.minmfixp[k], gdat.maxmfixp[k])
    
    if gdat.pntstype == 'lght':
        gdat.exprfwhm = 2. * retr_psfnwdth(gdat, gdat.exprpsfn, 0.5)
        gdat.stdvspatprio = amax(gdat.exprfwhm)
    if gdat.pntstype == 'lens':
        gdat.stdvspatprio = amax(gdat.exprpsfp)
    
    # proposals
    # parameters not subject to proposals
    gdat.indxfixpiact = setdiff1d(gdat.indxfixp, gdat.indxfixpprop)
    gdat.numbfixpiact = gdat.indxfixpiact.size
    gdat.indxiact = arange(gdat.numbfixpiact)
    
    gdat.lablproptype = array([])
    gdat.legdproptype = array([])
    gdat.nameproptype = array([])
   
    if gdat.probtran == None:
        if gdat.numbtrap > 0:
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
       
    cntr = tdpy.util.cntr()
    if gdat.numbtrap > 0.:
        
        if gdat.propwithsing:
            for k in gdat.indxstdp:    
                gdat.indxproptypewith = cntr.incr()
                gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{W_{%d}}$' % k)
                gdat.legdproptype = append(gdat.legdproptype, 'Within-model')
                gdat.nameproptype = append(gdat.nameproptype, gdat.namestdp[k])
        else:    
            gdat.indxproptypewith = cntr.incr()
            gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{W}$')
            gdat.legdproptype = append(gdat.legdproptype, 'Within-model')
            gdat.nameproptype = append(gdat.nameproptype, 'with')
    
        gdat.indxproptypebrth = cntr.incr()
        gdat.indxproptypedeth = cntr.incr()
        if gdat.probtran > 0.:
            # birth
            gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{B}$')
            gdat.legdproptype = append(gdat.legdproptype, 'Birth')
            gdat.nameproptype = append(gdat.nameproptype, 'brth')
            
            # death
            gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{D}$')
            gdat.legdproptype = append(gdat.legdproptype, 'Death')
            gdat.nameproptype = append(gdat.nameproptype, 'deth')
            
            if gdat.probbrde < 1.:
                # split
                gdat.indxproptypesplt = cntr.incr()
                gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{S}$')
                gdat.legdproptype = append(gdat.legdproptype, 'Split')
                gdat.nameproptype = append(gdat.nameproptype, 'splt')
                
                # merge
                gdat.indxproptypemerg = cntr.incr()
                gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{M}$')
                gdat.legdproptype = append(gdat.legdproptype, 'Merge')
                gdat.nameproptype = append(gdat.nameproptype, 'merg')
    
    gdat.numbproptype = gdat.nameproptype.size
    gdat.indxproptype = arange(gdat.numbproptype)

   
def setpfinl(gdat, boolinitsetp=False):

    # get count data
    if gdat.pixltype == 'cart':
        # temp
        gdat.indxxaximaxm, gdat.indxyaximaxm = tdpy.util.retr_indximagmaxm(gdat.datacnts[0, :, 0].reshape((gdat.numbsidecart, gdat.numbsidecart)))

    gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix
    if gdat.enerdiff:
        gdat.datafluxmean /= gdat.deltener

    # sanity checks
    # temp
    if (fabs(gdat.datacnts - rint(gdat.datacnts)) > 1e-3).any() and boolinitsetp:
        print 'Fractional counts!'

    if amin(gdat.datacnts) < 0. and boolinitsetp:
        print 'Negative counts!'


def retr_datatick(gdat):

    # data count limits
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
    

def retr_fromgdat(gdat, gdatmodi, strg, strgvarb, stdv=False, errr=False):
    
    if strgvarb == 'datacnts':
        varb = gdat.datacnts
    elif strg == 'this':
        varb = getattr(gdatmodi, strg + strgvarb)
    elif strg == 'post':
        if stdv:
            varb = getattr(gdat, 'stdv' + strgvarb)
        elif errr:
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
    oaxitype = getattr(gdat, strgpara + 'oaxitype')
    psfntype = getattr(gdat, strgpara + 'psfntype')
    maxmnumbpnts = getattr(gdat, strgpara + 'maxmnumbpnts') 
    numbback = len(getattr(gdat, strgpara + 'back'))
    specback = getattr(gdat, strgpara + 'specback')
    indxback = arange(numbback)

    # background
    ## number of background parameters
    numbbacp = 0
    for c in indxback:
        if specback[c] != None:
            numbbacp += 1
        else:
            numbbacp += gdat.numbener
   
    ## background parameter indices
    indxbacpback = []
    indxbackbacp = zeros(numbbacp, dtype=int)
    indxenerbacp = zeros(numbbacp, dtype=int)
    cntr = 0
    for c in indxback: 
        
        if specback[c] != None:
            indxbacpback.append(cntr)
            indxbackbacp[cntr] = c
            cntr += 1
        else:
            indxbacpback.append(cntr + gdat.indxener)
            for i in gdat.indxener:
                indxenerbacp[cntr] = i
                indxbackbacp[cntr] = c
                cntr += 1
   
    # total maximum number of elements
    maxmnumbpntstotl = sum(maxmnumbpnts)
    indxpntstotl = arange(maxmnumbpntstotl)
    maxmnumbpntscumr = cumsum(maxmnumbpnts)
    maxmnumbpntscuml = concatenate((array([0]), maxmnumbpntscumr[:-1]))
   
    # population index vector
    numbpopl = getattr(gdat, strgpara + 'numbpopl')
    indxpopl = arange(numbpopl, dtype=int) 

    liststrgcomp = [[] for l in indxpopl]
    listscalcomp = [[] for l in indxpopl]
    liststrgfeat = ['gang', 'aang', 'lgal', 'bgal', 'flux', 'spec', 'deltllik']
    
    liststrgfeatsign = ['flux', 'deltllik']

    if gdat.pntstype == 'lght':
        liststrgfeat += ['cnts']
    if gdat.pntstype == 'lens':
        liststrgfeat += ['distsour', 'dotpsour']
        liststrgfeatsign += ['distsour', 'dotpsour']

    liststrgfeatprio = [[] for l in indxpopl]
    if gdat.numbener > 1 and gdat.pntstype == 'lght':
        liststrgfeat += ['sind']
    for l in indxpopl:
        liststrgcomp[l] = ['lgal', 'bgal', 'flux']
        listscalcomp[l] = ['self', 'self', 'powr']
        if gdat.numbener > 1 and gdat.pntstype == 'lght':
            liststrgcomp[l] += ['sind']
            listscalcomp[l] += ['gaus']
        if spatdisttype[l] == 'gang':
            liststrgfeatprio[l] += ['gang', 'aang']
        if spatdisttype[l] == 'gaus':
            liststrgfeatprio[l] += ['lgal', 'bgal']
        elif spatdisttype[l] == 'disc' or spatdisttype[l] == 'unif':
            liststrgfeatprio[l] += ['lgal', 'bgal']
        liststrgfeatprio[l] += ['flux']
        if gdat.numbener > 1 and gdat.pntstype == 'lght':
            liststrgfeatprio[l] += ['sind']
            liststrgfeatsign += ['sind']
            if spectype[l] == 'curv':
                liststrgcomp[l] += ['curv']
                listscalcomp[l] += ['gaus']
                liststrgfeatprio[l] += ['curv']
                if not 'curv' in liststrgfeat:
                    liststrgfeat += ['curv']
            if spectype[l] == 'expo':
                liststrgcomp[l] += ['expo']
                listscalcomp[l] += ['lgau']
                liststrgfeatprio[l] += ['expo']
                if not 'expo' in liststrgfeat:
                    liststrgfeat += ['expo']
    
    # list of element parameters of all populations
    liststrgcomptotl = []
    for listsubb in liststrgcomp:
        for strg in listsubb:
            if not strg in liststrgcomptotl:
                liststrgcomptotl.append(strg)
    
    # list of element features that are not parameters
    liststrgfeatdiff = []
    for strgfeat in liststrgfeat:
        if not strgfeat in liststrgcomptotl:
            liststrgfeatdiff.append(strgfeat)

    liststrgfeatdefa = deepcopy(liststrgfeat)
    if not 'sind' in liststrgfeatdefa:
        liststrgfeatdefa += ['sind']
    if not 'curv' in liststrgfeatdefa:
        liststrgfeatdefa += ['curv']
    if not 'expo' in liststrgfeatdefa:
        liststrgfeatdefa += ['expo']

    cntr = tdpy.util.cntr()
    
    dicttemp = {}
    if maxmnumbpntstotl > 0:
        for l in indxpopl:
            dicttemp['indxfixpnumbpntspop%d' % l] = cntr.incr()
        for l in indxpopl:
            dicttemp['indxfixpmeanpntspop%d' % l] = cntr.incr()
    
        liststrgvarb = ['gangdistscal', 'bgaldistscal', 'spatdistcons', 'fluxdistslop', 'sinddistmean', 'sinddiststdv', \
                                                                                                    'curvdistmean', 'curvdiststdv', 'expodistmean', 'expodiststdv']
    
        # temp
        for strgvarb in liststrgvarb:
            strgtemp = 'indxfixp' + strgvarb
            dicttemp[strgtemp] = zeros(numbpopl, dtype=int) - 1
        dicttemp['indxfixpfluxdistnorm'] = zeros((numbpopl, gdat.numbfluxdistnorm), dtype=int)

        for l in range(numbpopl):
            for strg in liststrgfeatprio[l]:
                
                if strg == 'gang' and spatdisttype[l] == 'gang':
                    dicttemp['indxfixpgangdistscalpop%d' % l] = cntr.incr()
                    dicttemp['indxfixpgangdistscal'][l] = dicttemp['indxfixpgangdistscalpop%d' % l]
                if strg == 'bgal' and spatdisttype[l] == 'disc':
                    dicttemp['indxfixpbgaldistscalpop%d' % l] = cntr.incr()
                    dicttemp['indxfixpbgaldistscal'][l] = dicttemp['indxfixpbgaldistscalpop%d' % l]
                if strg == 'lgal' and spatdisttype[l] == 'gaus':
                    dicttemp['indxfixpspatdistconspop%d' % l] = cntr.incr()
                    dicttemp['indxfixpspatdistcons'][l] = dicttemp['indxfixpspatdistconspop%d' % l]
                if strg == 'flux':
                    if fluxdisttype[l] == 'powr':
                        dicttemp['indxfixpfluxdistsloppop%d' % l] = cntr.incr()
                        dicttemp['indxfixpfluxdistslop'][l] = dicttemp['indxfixpfluxdistsloppop%d' % l]
                    if fluxdisttype[l] == 'bind':
                        for k in range(gdat.numbfluxdistnorm):
                            dicttemp['indxfixpfluxdistnormbin%dpop%d' % (k, l)] = cntr.incr()
                            dicttemp['indxfixpfluxdistnorm'][l, k] = dicttemp['indxfixpfluxdistnormbin%dpop%d' % (k, l)]
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
        
    dicttemp['indxfixphypr'] = array(dicttemp['indxfixpdist'] +  dicttemp['indxfixpmeanpnts'])
   
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if psfntype == 'singgaus':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'singking':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'doubgaus':
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'gausking':
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'doubking':
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamtene%devt%d' % (i, m)] = cntr.incr()
            if oaxitype:
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
    
    if oaxitype:
        numbpsfpoaxi = 2
    else:
        numbpsfpoaxi = 0

    numbpsfptotl = numbpsfpform + numbpsfpoaxi
    
    if oaxitype:
        indxpsfponor = numbpsfpform + numbpsfptotl * arange(gdat.numbener * gdat.numbevtt)
        indxpsfpoind = numbpsfpform + numbpsfptotl * arange(gdat.numbener * gdat.numbevtt) + 1
    else:
        indxpsfponor = []
        indxpsfpoind = []

    numbpsfptotlevtt = gdat.numbevtt * numbpsfptotl
    numbpsfptotlener = gdat.numbener * numbpsfptotl
    numbpsfp = numbpsfptotl * gdat.numbener * gdat.numbevtt
    indxpsfpoaxi = arange(numbpsfpoaxi) 
    indxpsfpform = arange(numbpsfpform)
    indxpsfptotl = arange(numbpsfptotl)
   
    indxpsfp = arange(numbpsfp)
    indxfixppsfpinit = sort(dicttemp['indxfixppsfp'])[0]
    
    if oaxitype:
        indxfixppsfponor = indxfixppsfpinit + indxpsfponor
        indxfixppsfpoind = indxfixppsfpinit + indxpsfpoind
        indxfixppsfpoaxi = sort(concatenate((indxfixppsfponor, indxfixppsfpoind)))

    dicttemp['indxfixpbacp'] = []
    for c in indxback:
        if specback[c] != None:
            indx = cntr.incr()
            dicttemp['indxfixpbacpbac%d' % c] = indx
            dicttemp['indxfixpbacp'].append(indx)
        else:
            for i in gdat.indxener:
                indx = cntr.incr()
                dicttemp['indxfixpbacpbac%dene%d' % (c, i)] = indx
                dicttemp['indxfixpbacp'].append(indx)
    dicttemp['indxfixpbacp'] = array(dicttemp['indxfixpbacp'])
    dicttemp['indxfixphost'] = []
    dicttemp['indxfixpsour'] = []
    dicttemp['indxfixplenp'] = []
    dicttemp['indxfixpanglsour'] = []
    dicttemp['indxfixpanglhost'] = []
    dicttemp['indxfixpangllens'] = []
    dicttemp['indxfixpspecsour'] = []
    dicttemp['indxfixpspechost'] = []

    if gdat.pntstype == 'lens':
        dicttemp['indxfixplgalsour'] = cntr.incr()
        dicttemp['indxfixpbgalsour'] = cntr.incr()
        for i in gdat.indxener:
            dicttemp['indxfixpspecsourene%d' % i] = cntr.incr()
        dicttemp['indxfixpsizesour'] = cntr.incr()
        dicttemp['indxfixpellpsour'] = cntr.incr()
        dicttemp['indxfixpanglsour'] = cntr.incr()
        dicttemp['indxfixplgalhost'] = cntr.incr()
        dicttemp['indxfixpbgalhost'] = cntr.incr()
        for i in gdat.indxener:
            dicttemp['indxfixpspechostene%d' % i] = cntr.incr()
        dicttemp['indxfixpsizehost'] = cntr.incr()
        dicttemp['indxfixpbeinhost'] = cntr.incr()
        dicttemp['indxfixpellphost'] = cntr.incr()
        dicttemp['indxfixpanglhost'] = cntr.incr()
        dicttemp['indxfixpsherhost'] = cntr.incr()
        dicttemp['indxfixpsanghost'] = cntr.incr()
        dicttemp['indxfixpsour'] = []

        liststrgsour = ['lgalsour', 'bgalsour', 'specsour', 'sizesour', 'ellpsour', 'anglsour']
        liststrghost = ['lgalhost', 'bgalhost', 'spechost', 'sizehost', 'beinhost', 'ellphost', 'anglhost', 'sherhost', 'sanghost']
        for strg, valu in dicttemp.iteritems():
            
            if strg == 'indxfixpfluxdistnorm' or isinstance(valu, list) or isinstance(valu, ndarray):
                continue

            for strgtemp in liststrgsour:
                if strg[8:].startswith(strgtemp):
                    if isinstance(valu, list):
                        for valutemp in valu:
                            dicttemp['indxfixpsour'].append(valutemp)
                    else:
                        dicttemp['indxfixpsour'].append(valu)
            
            for strgtemp in liststrghost:
                if strg[8:].startswith(strgtemp):
                    if isinstance(valu, list):
                        for valutemp in valu:
                            dicttemp['indxfixphost'].append(valutemp)
                    else:
                        dicttemp['indxfixphost'].append(valu)
                
            if strg[8:].startswith('specsourene'):
                dicttemp['indxfixpspecsour'].append(valu)

            if strg[8:].startswith('spechostene'):
                dicttemp['indxfixpspechost'].append(valu)
            
            if valu in dicttemp['indxfixpsour'] or valu in dicttemp['indxfixphost']:
                dicttemp['indxfixplenp'].append(valu)
    
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

    # number of element parameters
    numbcomp = 3 + zeros(numbpopl, dtype=int)
    if gdat.numbener > 1:
        numbcomp += numbspep
    maxmnumbcomp = amax(numbcomp)

    indxcomp = []
    for l in indxpopl:
        indxcomp.append(arange(numbcomp[l]))

    # number of transdimensional parameters
    numbtrappopl = maxmnumbpnts * numbcomp
    numbtrapcumr = cumsum(numbtrappopl)
    numbtrapcuml = concatenate((array([0]), numbtrapcumr[:-1]))
    numbtrap = sum(numbtrappopl)
    
    # total number of parameters
    numbpara = numbfixp + numbtrap
    indxsamptrap = arange(numbfixp, numbpara)
    indxsampcompinit = indxsamptrap[0]
    indxpara = arange(numbpara)

    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr':
            setattr(gdat, strgpara + strg, valu)

    for strg, valu in dicttemp.iteritems():
        # temp
        if isinstance(valu, ndarray) and strg[12:16] != 'dist':
            valu = sort(valu)
        setattr(gdat, strgpara + strg, valu)


def setp_fixp(gdat, strgpara=''):
    
    listscalcomp = getattr(gdat, strgpara + 'listscalcomp')
    liststrgcomp = getattr(gdat, strgpara + 'liststrgcomp')
    numbcomp = getattr(gdat, strgpara + 'numbcomp')
    numbtrapcuml = getattr(gdat, strgpara + 'numbtrapcuml')
    numbtrapcumr = getattr(gdat, strgpara + 'numbtrapcumr')
    indxfixp = getattr(gdat, strgpara + 'indxfixp')
    numbfixp = getattr(gdat, strgpara + 'numbfixp')
    numbback = getattr(gdat, strgpara + 'numbback')
    numbtrap = getattr(gdat, strgpara + 'numbtrap')
    numbpara = getattr(gdat, strgpara + 'numbpara')
    numbpopl = getattr(gdat, strgpara + 'numbpopl')
    psfntype = getattr(gdat, strgpara + 'psfntype')
    numbpsfptotl = getattr(gdat, strgpara + 'numbpsfptotl')
    specback = getattr(gdat, strgpara + 'specback')
    indxbackbacp = getattr(gdat, strgpara + 'indxbackbacp')
    indxenerbacp = getattr(gdat, strgpara + 'indxenerbacp')

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
    
    for strg, k in gdat.__dict__.iteritems():
        
        if strgpara == 'true':
            strgvarb = strg[4:]
        else:
            strgvarb = strg

        if not isinstance(k, int) or not strgvarb.startswith('indxfixp') or strgvarb == 'indxfixppsfpinit':
            continue

        strgvarb = strgvarb[8:]
        
        namefixp[k] = strgvarb

        if strg[:-1].endswith('pop'):
            
            if numbpopl == 1:
                strgpopl = ''
                strgpoplcomm = ''
                strgpoplbind = strg[-5]
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
            
            if namefixp[k].startswith('spatdistconsp'):
                strgfixp[k] = r'$\gamma_{c%s}$' % strgpoplcomm
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('bgaldistscalp'):
                strgfixp[k] = r'$\gamma_{b%s}$' % strgpoplcomm
                scalfixp[k] = 'self'
            
            if namefixp[k].startswith('fluxdistslopp'):
                strgfixp[k] = r'$\alpha_{%s}$' % strgpopl
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('fluxdistnormbin'):
                strgfixp[k] = r'$p_{f,%s}$' % strgpoplbind
                scalfixp[k] = 'logt'
            
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
                n = k - getattr(gdat, strgpara + 'indxfixpsigcene0evt0')
                meanfixp[k] = gdat.exprpsfp[n]
                stdvfixp[k] = 0.1 * gdat.exprpsfp[n]
            else:
                if strgvarb.startswith('sig'):
                    scalfixp[k] = 'logt'
                if strgvarb.startswith('gam'):
                    scalfixp[k] = 'atan'
                if strgvarb.startswith('psff'):
                    scalfixp[k] = 'atan'
                if strgvarb.startswith('onor'):
                    scalfixp[k] = 'logt'
                if strgvarb.startswith('oind'):
                    scalfixp[k] = 'atan'
            
            # strings for PSF parameters
            if strgvarb.startswith('sig'):
                strgvarbtemp = '\sigma'
                factfixpplot[k] = gdat.anglfact
            if strgvarb.startswith('gam'):
                strgvarbtemp = '\gamma'
            if strgvarb.startswith('psff'):
                strgvarbtemp = 'f'
            if strgvarb.startswith('onor'):
                strgvarbtemp = 'a'
            if strgvarb.startswith('oind'):
                strgvarbtemp = 'b'
            if strgvarb.startswith('sig') and psfntype == 'doubgaus' or psfntype == 'gausking' or psfntype == 'doubking':
                strgcomptemp = 'c'
            elif strgvarb.startswith('gam') and psfntype == 'gausking' or psfntype == 'doubking':
                strgcomptemp = 't'
            else:
                strgcomptemp = ''
            if gdat.numbener > 1:
                indxenertemp = gdat.indxenerincl[((k - getattr(gdat, strgpara + 'indxfixpsigcene0evt0')) % (gdat.numbener * numbpsfptotl)) // numbpsfptotl]
                strgenertemp = '%s' % indxenertemp
            else:
                strgenertemp = ''
            if gdat.numbevtt > 1:
                indxevtttemp = gdat.indxevttincl[(k - getattr(gdat, strgpara + 'indxfixpsigcene0evt0')) // (gdat.numbener * numbpsfptotl)]
                strgevtttemp = '%s' % indxevtttemp
            else:
                strgevtttemp = ''
            strgfixp[k] = r'$%s^{%s}_{%s%s}$' % (strgvarbtemp, strgcomptemp, strgenertemp, strgevtttemp)
        
        if strgvarb.startswith('bacp'):
            
            try:
                indxfixpbacpinit
            except:
                indxfixpbacpinit = k
            
            c = indxbackbacp[k-indxfixpbacpinit]

            name = 'bacpbac%d' % c
            
            if specback[c] != None:
                strgenertemp = ''
            else:
                i = indxenerbacp[k-indxfixpbacpinit]
                if gdat.numbener > 1:
                    strgenertemp = '%d' % i
                else:
                    strgenertemp = ''
                name += 'ene%d' % i
            
            if numbback > 1:
                strgbacktemp = '%d' % c
            else:
                strgbacktemp = ''
            
            namefixp[k] = name
            strgfixp[k] = '$A_{%s%s}$' % (strgbacktemp, strgenertemp)
            scalfixp[k] = 'logt'
        
        if gdat.pntstype == 'lens':
            if k in getattr(gdat, strgpara + 'indxfixplenp'):
                if strgvarb == 'lgalsour':
                    strgfixp[k] = '$l_s$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strgvarb == 'bgalsour':
                    strgfixp[k] = '$b_s$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strgvarb.startswith('specsour'):
                    if gdat.numbener > 1:
                        strgfixp[k] = '$f_{s,%s}$' % strg[-1]
                    else:
                        strgfixp[k] = '$f_s$'
                    scalfixp[k] = 'logt'
                if strgvarb == 'sizesour':
                    strgfixp[k] = '$a_s$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                if strgvarb == 'ellpsour':
                    strgfixp[k] = r'$\epsilon_s$'
                    scalfixp[k] = 'self'
                if strgvarb == 'anglsour':
                    strgfixp[k] = r'$\phi_s$'
                    scalfixp[k] = 'self'
                if strgvarb == 'lgalhost':
                    strgfixp[k] = '$l_h$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strgvarb == 'bgalhost':
                    strgfixp[k] = '$b_h$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                if strgvarb.startswith('spechost'):
                    if gdat.numbener > 1:
                        strgfixp[k] = '$f_{h,%s}$' % strg[-1]
                    else:
                        strgfixp[k] = '$f_h$'
                    scalfixp[k] = 'logt'
                if strgvarb == 'sizehost':
                    strgfixp[k] = '$a_h$'
                    scalfixp[k] = 'logt'
                if strgvarb == 'beinhost':
                    strgfixp[k] = r'$\theta_{E,h}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                if strgvarb == 'ellphost':
                    strgfixp[k] = r'$\epsilon_h$'
                    scalfixp[k] = 'self'
                if strgvarb == 'anglhost':
                    strgfixp[k] = r'$\phi_h$'
                    scalfixp[k] = 'self'
                if strgvarb == 'sherhost':
                    strgfixp[k] = r'$\gamma_e$'
                    scalfixp[k] = 'self'
                if strgvarb == 'sanghost':
                    strgfixp[k] = r'$\phi_{\gamma}$'
                    scalfixp[k] = 'self'
        
        if scalfixp[k] == 'pois' or scalfixp[k] == 'self' or scalfixp[k] == 'logt' or scalfixp[k] == 'atan':
            if strgvarb.startswith('numbpnts'):
                minmfixp[k] = getattr(gdat, strgpara + 'minm' + strgvarb[:-4])[k]
                maxmfixp[k] = getattr(gdat, strgpara + 'maxm' + strgvarb[:-4])[k]
            else:
                minmfixp[k] = getattr(gdat, strgpara + 'minm' + strgvarb)
                maxmfixp[k] = getattr(gdat, strgpara + 'maxm' + strgvarb)
        if scalfixp[k] == 'gaus' or scalfixp[k] == 'eerr':
            meanfixp[k] = getattr(gdat, strgpara + 'mean' + strgvarb)
            stdvfixp[k] = getattr(gdat, strgpara + 'stdv' + strgvarb)
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
        
        if strgvarb.startswith('sig'):
            strgfixpunit[k] = strgfixp[k] + ' [%s]' % getattr(gdat, 'lablgangunit')
        else:
            strgfixpunit[k] = strgfixp[k]
    
    namepara = zeros(numbpara, dtype=object)
    scalpara = zeros(numbpara, dtype=object)
    namepara[indxfixp] = namefixp
    scalpara[indxfixp] = scalfixp
    for k in range(numbtrap):
        indxpopltemp = argmin(where(k // numbtrapcumr == 0))
        indxcomptemp = (k - numbtrapcuml[indxpopltemp]) % numbcomp[indxpopltemp]
        namepara[numbfixp+k] = liststrgcomp[indxpopltemp][indxcomptemp]
        scalpara[numbfixp+k] = listscalcomp[indxpopltemp][indxcomptemp]
    
    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr':
            setattr(gdat, strgpara + strg, valu)


def setp_truetemp(gdat, strgvarb, popl, ener, evtt, back):
    
    liststrgvarb = []
    if popl:
        for l in gdat.trueindxpopl:
            liststrgvarb.append(strgvarb + 'pop%d' % l)
    elif ener and evtt:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                liststrgvarb.append(strgvarb + 'ene%devt%d' % (i, m))
    elif ener and back:
        for c in gdat.trueindxback:
            liststrgvarb.append(strgvarb + 'bac%d' % c)
            for i in gdat.indxener:
                liststrgvarb.append(strgvarb + 'bac%dene%d' % (c, i))
    elif ener:
        for i in gdat.indxener:
            liststrgvarb.append(strgvarb + 'ene%d' % i)
    else:   
        liststrgvarb.append(strgvarb)

    return liststrgvarb


def setp_true(gdat, strgvarb, valu, popl=False, ener=False, evtt=False, back=False):
    
    liststrgvarb = setp_truetemp(gdat, strgvarb, popl, ener, evtt, back)

    for strgvarb in liststrgvarb:
        try:
            getattr(gdat, 'true' + strgvarb)
        except:
            setattr(gdat, 'true' + strgvarb, valu)
 

def setp_truedefa(gdat, strgvarb, listvalu, typelimt='minmmaxm', popl=False, ener=False, evtt=False, back=False):
    
    liststrgvarb = setp_truetemp(gdat, strgvarb, popl, ener, evtt, back)

    for strgvarb in liststrgvarb:

        if typelimt == 'minmmaxm':
            try:
                minmpara = getattr(gdat, 'trueminm' + strgvarb)
            except:
                setattr(gdat, 'trueminm' + strgvarb, listvalu[0])
            try:
                maxmpara = getattr(gdat, 'truemaxm' + strgvarb)
            except:
                setattr(gdat, 'truemaxm' + strgvarb, listvalu[1])
        else:
            try:
                meanpara = getattr(gdat, 'truemean' + strgvarb)
            except:
                setattr(gdat, 'truemean' + strgvarb, listvalu[0])
            try:
                stdvpara = getattr(gdat, 'truestdv' + strgvarb)
            except:
                setattr(gdat, 'truestdv' + strgvarb, listvalu[1])
 

def intp_sinc(gdat, lgal, bgal):

    intpsinc = 4. * gdat.numbsidepsfn**2 * sum(gdat.temppsfn * sinc(gdat.numbsidepsfn * (gdat.gridpsfnlgal + lgal) - gdat.gridpsfnlgal) * \
                                                                                                      sinc(gdat.numbsidepsfn * (gdat.gridpsfnbgal + bgal) - gdat.gridpsfnbgal))

    return intpsinc


def retr_fluxbrgt(gdat, lgal, bgal, flux):

    if lgal.size == 0:
        fluxbrgt = array([0.])
        fluxbrgtassc = array([0.])
    else:
        indxbrgt = argmax(flux)
        fluxbrgt = flux[indxbrgt]
        dist = retr_angldist(gdat, lgal, bgal, lgal[indxbrgt], bgal[indxbrgt])
        indxbrgtassc = where(dist < gdat.anglassc)[0]
        fluxbrgtassc = flux[indxbrgtassc]
        fluxbrgt = repeat(fluxbrgt, fluxbrgtassc.size)

    return fluxbrgt, fluxbrgtassc


def retr_indxoaxipnts(gdat, lgal, bgal):

    oaxi = retr_angldist(gdat, lgal, bgal, 0., 0.)
    indxoaxipnts = digitize(oaxi[0], gdat.binsoaxiopen) - 1
   
    return indxoaxipnts


def init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=None, indxevttplot=None, indxpoplplot=-1):

    if strg == 'this' or gdatmodi != None:
        pathfold = gdat.pathfram
    elif strg == 'true' or strg == '':
        pathfold = gdat.pathinit
    elif strg == 'post':
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
    titl = ''
    if indxenerplot != None and gdat.numbener > 1 and strgplot.endswith('cnts'):
        titl = gdat.strgener[indxenerplot]
    if indxevttplot != None and gdat.numbevtt > 1 and strgplot.endswith('cnts'):
        titl += ' ' + gdat.strgevtt[indxevttplot]
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


def retr_imag(gdat, axis, maps, strg, thisindxener=None, thisindxevtt=-1, cmap='Greys', vmin=None, vmax=None, scal=None, tdim=False):
    
    if scal == None:
        scal = gdat.scalmaps

    if vmin == None and vmax != None:
        vmin = -vmax
    
    draw_frambndr(gdat, axis)
    
    # flatten the array
    if tdim:
        if thisindxener == None:
            maps = maps.flatten()
        else:
            maps = maps.reshape((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    # take the relevant energy and PSF bins
    if thisindxener != None:
        if thisindxevtt == -1:
            maps = sum(maps[thisindxener, ...], axis=1)
        else:
            maps = maps[thisindxener, :, thisindxevtt]
    
    # project the map to 2D
    if gdat.pixltype == 'heal' and not tdim:
        maps = tdpy.util.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                            minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                            minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata)
    
    if gdat.pixltype == 'cart' or tdim:
        shap = [gdat.numbsidecart] + list(maps.shape)
        shap[1] = gdat.numbsidecart
        maps = maps.reshape(shap).swapaxes(0, 1)
   
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


def make_catllabl(gdat, strg, axis):
    
    if strg == 'post':
        colr = 'black'
        labl = 'Condensed Model %s' % gdat.strgelem
    else:
        colr = 'b'
        labl = 'Sample Model %s' % gdat.strgelem
    axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, label=labl, marker='+', linewidth=2, color=colr)
        
    if gdat.trueinfo:
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablhits, marker='x', linewidth=2, color='g')
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablbias, marker='*', linewidth=2, color='g', facecolor='none')
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, facecolor='none', \
                                                                                                    label=gdat.truelablmiss, marker='s', linewidth=2, color='g')
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
        
    axis.legend(bbox_to_anchor=[0.5, 1.15], loc='center', ncol=3)
        

def supr_fram(gdat, gdatmodi, strg, axis, indxpoplplot=-1):

    # true catalog
    if gdat.trueinfo:
        if indxpoplplot == -1:
            listindxpoplplot = gdat.trueindxpopl
        else:
            listindxpoplplot = [indxpoplplot]
        for l in listindxpoplplot:
        
            if gdat.numbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, gdat.truespec[l][0, gdat.indxenerfluxdist, :].flatten())
                lgal = copy(gdat.truelgal[l])
                bgal = copy(gdat.truebgal[l])
                numbpnts = int(gdat.truenumbpnts[l])
                
                if gdatmodi != None:
                   
                    ## associations
                    ### missed
                    indx = gdatmodi.trueindxpntsasscmiss[l]
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.truelablmiss, facecolor='none', \
                                                                                                                                marker='s', linewidth=2, color='g')
                    
                    ### biased
                    indx = gdatmodi.trueindxpntsasscbias[l]
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                                label=gdat.truelablbias, marker='*', linewidth=2, color='g', facecolor='none')
                    
                    ### hit
                    indx = gdatmodi.trueindxpntsasschits[l]
                    
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
    
    # model catalog
    if indxpoplplot == -1:
        listindxpoplplot = gdat.indxpopl
    else:
        listindxpoplplot = [indxpoplplot]
    for l in listindxpoplplot:
        if gdatmodi != None:
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
   
    # temp
    try:
        if strg == 'post':
            flux = array([gdat.dictglob['postelemdetr'][r]['flux'][0] for r in range(gdat.numbpntsdetr)])
            mrkrsize = retr_mrkrsize(gdat, flux)
            lgal = array([gdat.dictglob['postelemdetr'][r]['lgal'][0] for r in range(gdat.numbpntsdetr)])
            bgal = array([gdat.dictglob['postelemdetr'][r]['bgal'][0] for r in range(gdat.numbpntsdetr)])
            axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, label='Condensed', marker='+', linewidth=2, color='black')
            for r in len(gdat.numbpntsdetr):
                lgal = gdat.dictglob['listelemdetr'][r]['lgal']
                bgal = gdat.dictglob['listelemdetr'][r]['bgal']
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, marker='+', linewidth=2, color='black', alpha=0.3)
    except:
        pass



def retr_levi(listllik):
    
    minmlistllik = amin(listllik)
    levi = log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    return levi


def retr_infofromlevi(listllik, levi):
    
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


def retr_angldist(gdat, lgalfrst, bgalfrst, lgalseco, bgalseco):
    
    if gdat.pixltype == 'heal':
        dir1 = array([lgalfrst, bgalfrst])
        dir2 = array([lgalseco, bgalseco])
        angldist = hp.rotator.angdist(dir1, dir2)
    else:
        angldist = sqrt((lgalfrst - lgalseco)**2 + (bgalfrst - bgalseco)**2)

    return angldist


def retr_deflextr(gdat, sher, sang):
    
    factcosi = sher * cos(2. * sang)
    factsine = sher * cos(2. * sang)
    defllgal = factcosi * gdat.lgalgrid + factsine * gdat.bgalgrid
    deflbgal = factsine * gdat.lgalgrid - factcosi * gdat.bgalgrid
    
    return vstack((defllgal, deflbgal)).T 


def readfile(path):

    filepick = open(path + '.p', 'rb')
    filearry = h5py.File(path + '.h5', 'r')
    gdattemptemp = cPickle.load(filepick)
    
    for attr in filearry:
        setattr(gdattemptemp, attr, filearry[attr][()])

    filepick.close()
    filearry.close()
    
    return gdattemptemp


def readoutp(path):
    
    thisfile = h5py.File(path, 'r')
    gdattemp = tdpy.util.gdatstrt()
    for attr in thisfile:
        setattr(gdattemp, attr, thisfile[attr]) 
    thisfile.close()

    return gdattemp


def writfile(gdattemp, path):
    
    filepick = open(path + '.p', 'wb')
    filearry = h5py.File(path + '.h5', 'w')
    
    gdattemptemp = tdpy.util.gdatstrt()
    for attr, valu in gdattemp.__dict__.iteritems():
        if isinstance(valu, ndarray) and valu.dtype != dtype('O'):
            filearry.create_dataset(attr, data=valu)
        else:
            setattr(gdattemptemp, attr, valu)

    cPickle.dump(gdattemptemp, filepick, protocol=cPickle.HIGHEST_PROTOCOL)
    filepick.close()
    filearry.close()
   

def writoutp(gdat, path):

    thisfile = h5py.File(path, 'w')
    for attr, valu in gdat.__dict__.iteritems():
        if isinstance(valu, ndarray) and attr.startswith('list'):
            thisfile.create_dataset(attr, data=valu)
    thisfile.close()


def retr_deflcutf(angl, deflscal, anglscal, anglcutf=None, asym=False):

    fracscal = angl / anglscal
    
    fact = zeros_like(fracscal)
    indxlowr = where(fracscal < 1.)[0]
    indxuppr = where(fracscal > 1.)[0]
    fact[indxlowr] = arccosh(1. / fracscal[indxlowr]) / sqrt(1. - fracscal[indxlowr]**2)
    fact[indxuppr] = arccos(1. / fracscal[indxuppr]) / sqrt(fracscal[indxuppr]**2 - 1.)
    
    deflcutf = deflscal / fracscal / (1. - log(2.))
    if asym:
        deflcutf *= fact + log(fracscal / 2.)
    else:
        fraccutf = anglcutf / anglscal
        factcutf = fraccutf**2 / (fraccutf**2 + 1)**2 * ((fraccutf**2 + 1. + 2. * (fracscal**2 - 1.)) * fact + \
                pi * fraccutf + (fraccutf**2 - 1.) * log(fraccutf) + sqrt(fracscal**2 + fraccutf**2) * (-pi + (fraccutf**2 - 1.) / fraccutf * \
                log(fracscal / (sqrt(fracscal**2 + fraccutf**2) + fraccutf))))
        deflcutf *= factcutf
    
    return deflcutf


def retr_defl(gdat, lgal, bgal, bein, ellp, angl, rcor, indxpixltemp=None, cutf=False):
    
    if indxpixltemp == None:
        indxpixltemp = gdat.indxpixl
    
    # translate the grid
    lgaltran = gdat.lgalgrid[indxpixltemp] - lgal
    bgaltran = gdat.bgalgrid[indxpixltemp] - bgal
    
    if cutf:
        angl = sqrt(lgaltran**2 + bgaltran**2)
        defl = retr_deflcutf(angl, bein, anglscal=gdat.anglscal, anglcutf=gdat.anglcutf)
        
        defllgal = lgaltran / angl * defl
        deflbgal = bgaltran / angl * defl

    else:
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
    
    return vstack((defllgal, deflbgal)).T


def retr_lprifluxdist(gdat, gdatmodi, flux, sampvarb, l, minmflux, binsflux, fluxdisttype):
    
    if fluxdisttype[l] == 'powr':
        fluxdistslop = sampvarb[gdat.indxfixpfluxdistslop[l]]
        lprifluxdist = sum(log(pdfn_powr(flux, minmflux, gdat.maxmflux, fluxdistslop)))
    if fluxdisttype[l] == 'bind':
        fluxdistnorm = sampvarb[gdat.indxfixpfluxdistnorm[l, :]]
        lprifluxdist = sum(log(pdfn_bind(flux, minmflux, gdat.maxmflux, binsflux, fluxdistnorm)))

    return lprifluxdist


def retr_lprisinddist(gdat, gdatmodi, sind, sampvarb, l, sinddisttype):
    
    if sinddisttype[l] == 'gaus':
        lprisinddist = sum(log(pdfn_gaus(sind, sampvarb[gdat.indxfixpsinddistmean[l]], sampvarb[gdat.indxfixpsinddiststdv[l]]))) 
    if sinddisttype[l] == 'atan':
        lprisinddist = sum(log(pdfn_atan(sind, gdat.minmsind, gdat.maxmsind))) 
    
    return lprisinddist


def retr_lpricurvdist(gdat, gdatmodi, curv, sampvarb, l):
    
    lpri = sum(log(pdfn_gaus(curv, sampvarb[gdat.indxfixpcurvdistmean[l]], sampvarb[gdat.indxfixpcurvdiststdv[l]])))
    
    return lpri


def traptdim(gdat, arry):
    
    s1 = arry[0, 0] + arry[-1, 0] + arry[0, -1] + arry[-1, -1]
    s2 = sum(arry[1:-1, 0]) + sum(arry[1:-1, -1]) + sum(arry[0, 1:-1]) + sum(arry[-1, 1:-1])
    s3 = sum(arry[1:-1, 1:-1])
    summ = (s1 + 2*s2 + 4*s3) * gdat.apixmodl
    
    return summ


def retr_spatprio(gdat, spatdistcons, pdfnspatpriotemp):
    
    pdfnspatprio = spatdistcons + pdfnspatpriotemp
    summ = traptdim(gdat, pdfnspatprio)
    pdfnspatprio /= summ
    lpdfspatprio = log(pdfnspatprio)
    lpdfspatprioobjt = sp.interpolate.RectBivariateSpline(gdat.binsbgalcart, gdat.binslgalcart, lpdfspatprio)
    
    return lpdfspatprio, lpdfspatprioobjt


def retr_lpriexpodist(gdat, gdatmodi, expo, sampvarb, l):
    
    lpri = sum(log(pdfn_gaus(expo, sampvarb[gdat.indxfixpexpodistmean[l]], sampvarb[gdat.indxfixpexpodiststdv[l]])))
    
    return lpri


def proc_samp(gdat, gdatmodi, strg, raww=False, fast=False, lprionly=False):

    # initial
    if gdatmodi != None:
        timeinit = gdat.functime()

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
    else:
        strgtype = ''
    indxpopl = getattr(gdat, strgtype + 'indxpopl')
    spectype = getattr(gdat, strgtype + 'spectype')
    fluxdisttype = getattr(gdat, strgtype + 'fluxdisttype')
    if gdat.numbener > 1:
        sinddisttype = getattr(gdat, strgtype + 'sinddisttype')
    spatdisttype = getattr(gdat, strgtype + 'spatdisttype')
    if gdat.pntstype == 'lght':
        oaxitype = getattr(gdat, strgtype + 'oaxitype')

    # temp
    minmflux = getattr(gdat, strgtype + 'minmflux')
    binsflux = getattr(gdat, strgtype + 'binsflux')
    
    # common dictionary
    dicttemp = {}
           
    # grab the sample vector
    sampvarb = getattr(gdatobjt, strg + 'sampvarb')
    
    psfp = sampvarb[getattr(gdat, strgtype + 'indxfixppsfp')]
    bacp = sampvarb[getattr(gdat, strgtype + 'indxfixpbacp')]
    
    indxpntsfull = list(getattr(gdatobjt, strg + 'indxpntsfull'))
    dicttemp['indxsamplgal'], dicttemp['indxsampbgal'], dicttemp['indxsampflux'], dicttemp['indxsampsind'], dicttemp['indxsampcurv'], \
																				dicttemp['indxsampexpo'], dicttemp['indxsampcomp'] = retr_indxsampcomp(gdat, indxpntsfull, spectype)
    for strgcomp in gdat.liststrgcomptotl:
        setattr(gdatobjt, strg + 'indxsamp' + strgcomp, dicttemp['indxsamp' + strgcomp])
    setattr(gdatobjt, strg + 'indxsampcomp', dicttemp['indxsampcomp'])
    
    if strg == 'next' and gdat.verbtype > 1:
        show_samp(gdat, gdatmodi)
    
    numbpnts = sampvarb[gdat.indxfixpnumbpnts].astype(int)
    numbpopl = numbpnts.size
    
    for strgfeat in gdat.liststrgfeatdefa:
        dicttemp[strgfeat] = [[] for l in range(numbpopl)]
    
    for l in range(numbpopl):
    	for strgcomp in gdat.liststrgcomp[l]:
            dicttemp[strgcomp][l] = sampvarb[dicttemp['indxsamp' + strgcomp][l]]
    
    for l in range(numbpopl):
        dicttemp['spec'][l] = retr_spec(gdat, dicttemp['flux'][l], dicttemp['sind'][l], dicttemp['curv'][l], dicttemp['expo'][l], spectype=spectype[l])
   
    lgalconc = concatenate(dicttemp['lgal'])
    bgalconc = concatenate(dicttemp['bgal'])
    specconc = concatenate(dicttemp['spec'], axis=1)
    numbpntsconc = lgalconc.size

    if gdatmodi != None:
        timefinl = gdat.functime()
        gdatmodi.thischroproc[0] = timefinl - timeinit
    
    # log-prior
    if gdatmodi != None:
        timeinit = gdat.functime()

    if gdat.numbtrap > 0:
    
        if 'gaus' in spatdisttype:
            spatdistcons = sampvarb[getattr(gdat, strgtype + 'indxfixpspatdistcons')]
            pdfnspatpriotemp = getattr(gdat, strgtype + 'pdfnspatpriotemp')
            lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, spatdistcons, pdfnspatpriotemp)
            lpdfspatpriointp = lpdfspatprioobjt(gdat.bgalcart, gdat.lgalcart)
            
            setattr(gdatobjt, strg + 'lpdfspatpriointp', lpdfspatpriointp)
            setattr(gdatobjt, strg + 'lpdfspatprioobjt', lpdfspatprioobjt)
        
        meanpnts = sampvarb[gdat.indxfixpmeanpnts]
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
                lpri[1+2*gdat.numbpopl+l] = -numbpnts[l] * log(2. * pi) 
            elif spatdisttype[l] == 'gaus':
                lpri[1+gdat.numbpopl+l] = sum(lpdfspatprioobjt(dicttemp['bgal'][l], dicttemp['lgal'][l], grid=False)) 
            lpri[1+3*gdat.numbpopl+l] = retr_lprifluxdist(gdat, gdatmodi, dicttemp['flux'][l], sampvarb, l, minmflux, binsflux, fluxdisttype)
            if gdat.numbener > 1:
                lpri[1+4*gdat.numbpopl+l] = retr_lprisinddist(gdat, gdatmodi, dicttemp['sind'][l], sampvarb, l, sinddisttype)
                if gdat.spectype[l] == 'curv':
                    lpri[1+5*gdat.numbpopl+l] = retr_lpricurvdist(gdat, gdatmodi, dicttemp['curv'][l], sampvarb, l)
                if gdat.spectype[l] == 'expo':
                    lpri[1+6*gdat.numbpopl+l] = retr_lpriexpodist(gdat, gdatmodi, dicttemp['expo'][l], sampvarb, l)
       
        lpritotl = sum(lpri)
        
        if strg == 'next' and (gdatmodi.propbrth or gdatmodi.propdeth):
            
            gdatmodi.thislpau = zeros(gdat.numblpau)
            
            if gdatmodi.propbrth:
                sampvarbtemp = gdatmodi.nextsampvarb
            if gdatmodi.propdeth:
                sampvarbtemp = gdatmodi.thissampvarb
        
            for l in gdat.indxpopl:
                
                gdatmodi.thislpau[gdat.maxmnumbcomp*l+0] = -log(2. * gdat.maxmgang)
                gdatmodi.thislpau[gdat.maxmnumbcomp*l+1] = -log(2. * gdat.maxmgang)
                gdatmodi.thislpau[gdat.maxmnumbcomp*l+2] = retr_lprifluxdist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[2]], \
                                                                                           gdatmodi.thissampvarb, gdatmodi.indxpoplmodi, minmflux, binsflux, fluxdisttype)
                if gdat.numbener > 1:
                    gdatmodi.thislpau[gdat.maxmnumbcomp*l+3] = retr_lprisinddist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[3]], \
                                                                                                 gdatmodi.thissampvarb, gdatmodi.indxpoplmodi, sinddisttype)
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'curv':
                        gdatmodi.thislpau[gdat.maxmnumbcomp*l+4] = retr_lpricurvdist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[4]], \
                                                                                                                gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                    if gdat.spectype[gdatmodi.indxpoplmodi] == 'expo':
                        gdatmodi.thislpau[gdat.maxmnumbcomp*l+4] = retr_lpriexpodist(gdat, gdatmodi, sampvarbtemp[gdatmodi.indxsamptran[4]], \
                                                                                                            gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
            
            if gdatmodi.propbrth:
                gdatmodi.thislpau *= -1.
            
            gdatmodi.thislpautotl = sum(gdatmodi.thislpau)
    
    setattr(gdatobjt, strg + 'lpritotl', lpritotl) 
    setattr(gdatobjt, strg + 'lpri', lpri)
    
    if gdatmodi != None:
        timefinl = gdat.functime()
        gdatmodi.thischroproc[1] = timefinl - timeinit
    
    ### loglikelihood
    if gdatmodi != None:
        timeinit = gdat.functime()

    if not lprionly:
        
        # process a sample vector and the occupancy list to calculate secondary variables
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
           
            ## host halo deflection
            defl = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost, 0.)
            
            ## external shear
            deflextr = retr_deflextr(gdat, sherhost, sanghost)
            defl += deflextr
        
            ## PSF kernel
            psfnkern = []
            for i in gdat.indxener:
                psfnkern.append(AiryDisk2DKernel(psfp[i] / gdat.sizepixl))
            
        if gdat.pntstype == 'lght':
            ## PSF off-axis factor
            if oaxitype:
                onor = sampvarb[getattr(gdat, strgtype + 'indxfixppsfponor')]
                oind = sampvarb[getattr(gdat, strgtype + 'indxfixppsfpoind')]
                factoaxi = retr_factoaxi(gdat, gdat.binsoaxiplot, onor, oind)
        
            psfntype = getattr(gdat, strgtype + 'psfntype')
            
            psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsanglplot, psfntype, gdat.binsoaxiplot, oaxitype)
            
            if oaxitype:
                psfnintp = []
                for p in gdat.indxoaxi:
                    psfnintp.append(interp1d_pick(gdat.binsanglplot, psfn[:, :, :, p], axis=1))
            else:
                psfnintp = interp1d_pick(gdat.binsanglplot, psfn, axis=1)
            setattr(gdatobjt, strg + 'psfnintp', psfnintp)
        
        
        if gdat.pntstype == 'lens':
            
            ## subhalos
            if numbpntsconc > 0 and not raww:
                deflelem = zeros((gdat.numbpixl, 2))
                for k in range(numbpntsconc):
                    # temp -- fix truncation
                    if gdat.evalcirc == 'full':
                        indxpixltemp = gdat.indxpixl
                    else:
                        indxpixlpnts = retr_indxpixl(gdat, bgalconc[k], lgalconc[k])
                        indxproxtemp = digitize(specconc[0, k], gdat.binsprox) - 1
                        indxpixltemp = gdat.indxpixlprox[indxproxtemp][indxpixlpnts]
                        if isinstance(indxpixltemp, int):
                            indxpixltemp = gdat.indxpixl
                    deflelem[indxpixltemp, :] += retr_defl(gdat, lgalconc[k], bgalconc[k], specconc[0, k], 0., 0., 0., indxpixltemp=indxpixltemp, cutf=True)
                defl += deflelem
            defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
            
            # lensed image
            lensflux = retr_mapssers(gdat, gdat.lgalgridcart - defl[:, :, 0], gdat.bgalgridcart - defl[:, :, 1], lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
            # host emission
            hostfluxmaps = retr_mapssers(gdat, gdat.lgalgridcart, gdat.bgalgridcart, lgalhost, bgalhost, spechost, sizehost, ellphost, anglhost)
           
            # total emission
            modlfluxuncv = lensflux + bacp * gdat.backfluxcart[0] + hostfluxmaps
            
            # convolve the lensed image with the PSF
            # temp
            # modlflux = modlfluxuncv.reshape((gdat.numbener, gdat.numbpixl, 1))
            modlflux = empty_like(gdat.expo)
            for i in gdat.indxener:
                modlflux[i, :, 0] = convolve_fft(modlfluxuncv[i, :, :, 0], psfnkern[i]).flatten()
            setattr(gdatobjt, strg + 'defl', defl)
            
        if gdat.pntstype == 'lght':
            
            ### PS flux map
            pntsflux = retr_pntsflux(gdat, lgalconc, bgalconc, specconc, psfnintp, oaxitype, evalcirc=gdat.evalcirc)
        
            setattr(gdatobjt, strg + 'pntsflux', pntsflux)
            
            ### model flux map
            modlflux = retr_mapslght(gdat, bacp, pntsflux, gdat.indxcube)
        
        ### count map
        modlcnts = retr_cntsmaps(gdat, modlflux)
        setattr(gdatobjt, strg + 'modlcnts', modlcnts)
        
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
            else:
                retr_datatick(gdat)
        
        ### log-likelihood
        llik = retr_llik(gdat, modlcnts)
        lliktotl = sum(llik)
    
        setattr(gdatobjt, strg + 'llik', llik) 
    else:
        lliktotl = 0.

    setattr(gdatobjt, strg + 'lliktotl', lliktotl) 
        
    if gdatmodi != None:
        timefinl = gdat.functime()
        gdatmodi.thischroproc[2] = timefinl - timeinit
    
    if lprionly:
        return
    else:
        lpostotl = lpritotl + lliktotl
        setattr(gdatobjt, strg + 'lpostotl', lpostotl) 

    if strg == 'next':
        setattr(gdatmodi, 'thislpriprop', lpri)
    
    if fast:
        return

    ## tertiary variables that are not needed for likelihood evaluation
    if strg != 'next':
       
        setattr(gdatobjt, strg + 'psfp', psfp)
        resicnts = gdat.datacnts - modlcnts
        setattr(gdatobjt, strg + 'resicnts', resicnts)
        
        ## derived variables
        if gdat.pntstype == 'lght':
            
            specplot = [[] for l in gdat.trueindxpopl]
            if l in indxpopl:
                #gdat.truespecplot[l] = empty((gdat.numbenerplot, gdat.truenumbpnts[l]))
                specplot[l] = retr_spec(gdat, dicttemp['flux'][l], dicttemp['sind'][l], dicttemp['curv'][l], dicttemp['expo'][l], spectype[l], plot=True)
            setattr(gdatobjt, strg + 'specplot', specplot)
             
            if oaxitype:
                setattr(gdatobjt, strg + 'factoaxi', factoaxi)
           
            setattr(gdatobjt, strg + 'psfn', psfn)

            ### PSF FWHM
            fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            setattr(gdatobjt, strg + 'fwhm', fwhm)
            
	    	### mean PS flux map 
            pntsfluxmean = sum(sum(pntsflux * gdat.expo, 2), 1) / sum(sum(gdat.expo, 2), 1)
            setattr(gdatobjt, strg + 'pntsfluxmean', pntsfluxmean)

            if gdat.calcerrr and gdat.numbtrap > 0:
                pntsflux = retr_pntsflux(gdat, lgalconc, bgalconc, specconc, psfnintp, gdat.oaxitype, evalcirc=False)
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

            ### sort with respect to deflection at scale radius
            if numbpntsconc > 0:
                indxpntssortbrgt = argsort(dicttemp['flux'][0])[::-1]
                lgalsort = dicttemp['lgal'][0][indxpntssortbrgt][:numbpntsconc]
                bgalsort = dicttemp['bgal'][0][indxpntssortbrgt][:numbpntsconc]
                beinsort = dicttemp['flux'][0][indxpntssortbrgt][:numbpntsconc]
    
            ### mass budget
            masspnts = gdat.massfromdeflscal * dicttemp['flux'][0]
            masspntstotl = array([sum(masspnts)])
            masshost = array([gdat.massfrombein * beinhost**2])
            fracsubh = masspntstotl / masshost
            
            setattr(gdatobjt, strg + 'masspntstotl', masspntstotl)
            setattr(gdatobjt, strg + 'masshost', masshost)
            setattr(gdatobjt, strg + 'fracsubh', fracsubh)

            deflsing = zeros((gdat.numbpixl, 2, gdat.numbdeflsingplot))
            numbdeflsing = min(gdat.numbdeflpntsplot, numbpntsconc) + 2
            for k in range(numbdeflsing):
                if k == 0:
                    deflsing[:, :, k] = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost, 0.)
                elif k == 1:
                    deflsing[:, :, k] = deflextr
                else:
                    if gdat.evalcirc == 'full':
                        indxpixltemp = gdat.indxpixl
                    else:
                        indxpixlpnts = retr_indxpixl(gdat, bgalsort[k-2], lgalsort[k-2])
                        indxproxtemp = digitize(beinsort[k-2], gdat.binsprox) - 1
                        indxpixltemp = gdat.indxpixlprox[indxproxtemp][indxpixlpnts]
                        if isinstance(indxpixltemp, int):
                            indxpixltemp = gdat.indxpixl
                    deflsing[indxpixltemp, :, k] = retr_defl(gdat, lgalsort[k-2], bgalsort[k-2], beinsort[k-2], 0., 0., 0., indxpixltemp=indxpixltemp, cutf=True)
            deflsing = deflsing.reshape((gdat.numbsidecart, gdat.numbsidecart, 2, gdat.numbdeflsingplot))

            ### convergence
            deflelem = deflelem.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
        
            conv = retr_conv(gdat, defl) 
            convelem = retr_conv(gdat, deflelem) 
            
            convpsec = retr_psec(gdat, conv)
            convpsecelem = retr_psec(gdat, conv)
            convpsecodim = retr_psecodim(gdat, convpsec) 
            convpsecelemodim = retr_psecodim(gdat, convpsecelem) 
            
            histdefl = histogram(defl, bins=gdat.binsdeflplot)[0]
            setattr(gdatobjt, strg + 'conv', conv)
            setattr(gdatobjt, strg + 'convelem', convelem)
            setattr(gdatobjt, strg + 'convpsec', convpsec)
            setattr(gdatobjt, strg + 'convpsecelem', convpsecelem)
            setattr(gdatobjt, strg + 'convpsecodim', convpsecodim)
            setattr(gdatobjt, strg + 'convpsecelemodim', convpsecelemodim)
            
            setattr(gdatobjt, strg + 'histdefl', histdefl)
            setattr(gdatobjt, strg + 'deflsing', deflsing)
     
        ## element features
        if gdat.numbtrap > 0:
            
            ### derived quantities
            for l in range(numbpopl):
                #### radial and angular coordinates
                dicttemp['gang'][l] = retr_gang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                dicttemp['aang'][l] = retr_aang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                
                if gdat.pntstype == 'lght':
                    #### number of expected counts
                    dicttemp['cnts'][l] = retr_pntscnts(gdat, dicttemp['lgal'][l], dicttemp['bgal'][l], dicttemp['spec'][l])
                
            #### delta log-likelihood
            gdatmoditemp = tdpy.util.gdatstrt()
            for l in range(numbpopl):
                dicttemp['deltllik'][l] = zeros(numbpnts[l])
                
                for k in range(numbpnts[l]):
                    if gdat.pntstype == 'lens':
                        defltemp = copy(defl.reshape((gdat.numbpixl, 2)))
                        defltemp[indxpixltemp, :] -= retr_defl(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['flux'][l][k], \
                                                                                                                0., 0., 0., indxpixltemp=indxpixltemp, cutf=True)
                        defltemp = defltemp.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
                        lensflux = retr_mapssers(gdat, gdat.lgalgridcart - defltemp[:, :, 0], gdat.bgalgridcart - defltemp[:, :, 1], \
                                                                                                            lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
                        modlfluxuncv = lensflux + bacp * gdat.backfluxcart[0] + hostfluxmaps
                        modlflux = empty_like(gdat.expo)
                        for i in gdat.indxener:
                            modlflux[i, :, 0] = convolve_fft(modlfluxuncv[i, :, :, 0], psfnkern[i]).flatten()
                    if gdat.pntstype == 'lght':
                        pntsfluxtemp = copy(pntsflux)
                        pntsfluxtemp -= retr_pntsflux(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['spec'][l][:, k], psfnintp, oaxitype, evalcirc=gdat.evalcirc)
                        modlflux = retr_mapslght(gdat, bacp, pntsfluxtemp, gdat.indxcube)
                    modlcnts = retr_cntsmaps(gdat, modlflux)
                    nextllik = retr_llik(gdat, modlcnts)
                    nextlliktotl = sum(nextllik)
                    dicttemp['deltllik'][l][k] = arcsinh(lliktotl - nextlliktotl) / log(10.)
                 
            # temp
            if False and gdat.pntstype == 'lght':
                ### number of background counts per PSF
                cntsbackfwhm = retr_cntsbackfwhm(gdat, bacp, fwhm)
            
                ### number of counts and standard deviation of each PS
                sigm = []
                for l in gdat.indxpopl:
                    # temp -- zero exposure pixels will give zero counts
                    if gdat.oaxitype:
                        sigmtemp = retr_sigm(gdat, dicttemp['cnts'][l], cntsbackfwhm, lgal=dicttemp['lgal'][l], bgal=dicttemp['bgal'][l])
                    else:
                        sigmtemp = retr_sigm(gdat, dicttemp['cnts'][l], cntsbackfwhm)
                    sigm.append(sigmtemp)
                
            if gdat.pntstype == 'lens':
                #### distance to the source
                for l in range(numbpopl):
                    dicttemp['distsour'][l] = retr_angldist(gdat, dicttemp['lgal'][l],  dicttemp['bgal'][l], sampvarb[gdat.indxfixplgalsour], sampvarb[gdat.indxfixpbgalsour])
                
                lensfluxmean = mean(sum(lensflux, 3), 0)
                
                #### dot product with the source flux gradient
                for l in range(numbpopl):
                    grad = dstack((gradient(lensfluxmean, gdat.sizepixl, axis=0), gradient(lensfluxmean, gdat.sizepixl, axis=1)))
                    dicttemp['dotpsour'][l] = empty(numbpnts[l])
                    deflsing = zeros((gdat.numbpixl, 2, numbpnts[l]))
                    for k in range(numbpnts[l]):
                        deflsing[indxpixltemp, :, k] = retr_defl(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['flux'][l][k], 0., 0., 0., \
                                                                                                                                        indxpixltemp=indxpixltemp, cutf=True)
                        deflsingtemp = deflsing[:, :, k].reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
                        dotpsourtemp = sum(grad * deflsingtemp, 2)
                        #fact = sqrt(grad[:, :, 0]**2 + grad[:, :, 1]**2) * sqrt(deflsingtemp[:, :, 0]**2 + deflsingtemp[:, :, 1]**2)
                        #indxfact = where(fact > 0.)
                        #dicttemp['dotpsour'][l][k] = mean(dotpsourtemp[indxfact] / fact[indxfact])
                        dicttemp['dotpsour'][l][k] = sum(fabs(dotpsourtemp))
            
            ### distribution of element parameters and features
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
                
            ### priors on element parameters and features
            for strgfeat in gdat.liststrgfeat:
                dicttemp[strgfeat + 'histprio'] = empty((numbpopl, gdat.numbbinsplotprio))
            
                # temp -- this does not work for mismodeling, need strg +
                #minmplot = getattr(gdat, 'minm' + strgfeat + 'plot')
                #maxmplot = getattr(gdat, 'maxm' + strgfeat + 'plot')
                if strgfeat == 'flux':
                    minm = getattr(gdat, strgtype + 'minm' + strgfeat)
                    maxm = getattr(gdat, strgtype + 'maxm' + strgfeat)
                    bins = getattr(gdat, strgtype + 'bins' + strgfeat)
                
                for l in range(numbpopl):
                    if strgfeat in gdat.liststrgfeatprio[l]:
                        
                        xdat = getattr(gdat, 'mean' + strgfeat + 'plotprio')
                        deltprio = getattr(gdat, 'delt' + strgfeat + 'plotprio')
                        delt = getattr(gdat, 'delt' + strgfeat + 'plot')
                        
                        booltemp = False
                        if strgfeat == 'gang' and spatdisttype[l] == 'gang' or strgfeat == 'bgal' and spatdisttype[l] == 'disc': 
                            scal = sampvarb[getattr(gdat, strgtype + 'indxfixp' + strgfeat + 'distscal')[l]]
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
                            slop = sampvarb[getattr(gdat, strgtype + 'indxfixp' + strgfeat + 'distslop')[l]]
                            pdfn = pdfn_powr(xdat, minm, maxm, slop)
                            booltemp = True
                        elif strgfeat == 'flux' and fluxdisttype[l] == 'bind':
                            fluxdistnorm = sampvarb[gdat.indxfixpfluxdistnorm[l, :]]
                            pdfn = pdfn_bind(xdat, minm, maxm, bins, fluxdistnorm)
                            booltemp = True
                        elif strgfeat == 'sind' or strgfeat == 'curv' and spectype[l] == 'curv' or strgfeat == 'expo' and spectype[l] == 'expo':
                            # this does not work for mismodeling
                            meanvarb = sampvarb[getattr(gdat, strgtype + 'indxfixp' + strgfeat + 'distmean')[l]]
                            stdv = sampvarb[getattr(gdat, strgtype + 'indxfixp' + strgfeat + 'diststdv')[l]]
                            if strgfeat == 'expo' and spectype[l] == 'expo':
                                pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                            else:
                                pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                            booltemp = True
                        
                        if False and strgfeat == 'sind':
                            print 'pdfn'
                            print pdfn
                            print 'deltprio'
                            print deltprio
                            print 'delt'
                            print delt
                            print 'meanpnts'
                            print meanpnts
                            print 

                        if booltemp:
                            dicttemp[strgfeat + 'histprio'][l, :] = meanpnts[l] * pdfn * deltprio * delt[0] / deltprio[0]
        
        for strgfeat in gdat.liststrgfeatdiff:
            if strgfeat != 'spec':
                setattr(gdatobjt, strg + strgfeat, dicttemp[strgfeat])

        for strgfeat in gdat.liststrgfeat:
            setattr(gdatobjt, strg + strgfeat + 'hist', dicttemp[strgfeat + 'hist'])
            setattr(gdatobjt, strg + strgfeat + 'histprio', dicttemp[strgfeat + 'histprio'])

        ### PS indices to compare with the reference catalog
        if strg == 'this' and gdat.numbtrap > 0 and gdat.trueinfo:
            
            if gdat.pntstype == 'lens' and gdat.trueinfo and gdat.datatype == 'mock':
                gdatmodi.thisdeflsingresi = gdatmodi.thisdeflsing - gdat.truedeflsing
                gdatmodi.thisdeflresi = gdatmodi.thisdefl - gdat.truedefl
                gdatmodi.thisdeflcomp = 1. - sum(gdatmodi.thisdefl * gdat.truedefl, axis=2) / sqrt(sum(gdatmodi.thisdefl**2, axis=2)) / sqrt(sum(gdat.truedefl**2, axis=2))
           
        if gdatmodi != None:
           
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
        gdatmodi.trueindxpntsasscmiss = [[] for l in gdat.trueindxpopl]
        gdatmodi.trueindxpntsasscbias = [[] for l in gdat.trueindxpopl]
        gdatmodi.trueindxpntsasschits = [[] for l in gdat.trueindxpopl]
        gdatmodi.trueindxpntsasscmult = [[] for l in gdat.trueindxpopl]
        gdatmodi.thisspecassc = [[] for l in gdat.trueindxpopl] 
        for l in gdat.trueindxpopl:
    
            indxmodlpnts = zeros(gdat.truenumbpnts[l], dtype=int) - 1
            specassc = zeros((gdat.numbener, gdat.truenumbpnts[l]))
            numbassc = zeros(gdat.truenumbpnts[l])
            metrassc = zeros(gdat.truenumbpnts[l]) + 3 * gdat.maxmgang
        
            for k in range(gdatmodi.thissampvarb[gdat.indxfixpnumbpnts[l]].astype(int)):
               
                # determine which true PSs satisfy the match criterion
                if gdat.asscmetrtype == 'dist':
                    metr = retr_angldist(gdat, gdat.truelgal[l], gdat.truebgal[l], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k])
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
                    gdatmodi.trueindxpntsasscmiss[l].append(k)
                else:
                    if numbassc[k] > 1:
                        gdatmodi.trueindxpntsasscmult[l].append(k)
            
                    ## check whether the flux of the associated model point source matches well with the flux of the deterministic point source
                    boolbias = specassc[gdat.indxenerfluxdist[0], k] > fluxbias[1, gdat.indxenerfluxdist[0], k] or \
                                                                                            specassc[gdat.indxenerfluxdist[0], k] < fluxbias[0, gdat.indxenerfluxdist[0], k]
                    if boolbias:
                        gdatmodi.trueindxpntsasscbias[l].append(k)
                    else:
                        gdatmodi.trueindxpntsasschits[l].append(k)
            
            gdatmodi.thisspecassc[l] = zeros((gdat.numbener, gdat.truenumbpnts[l]))
            temp = where(indxmodlpnts >= 0)[0]
            gdatmodi.thisspecassc[l][:, temp] = dicttemp['spec'][l][:, indxmodlpnts[temp]]
    
    
def retr_info(pdfnpost, pdfnprio):
    
    info = pdfnpost * log(pdfnpost / pdfnprio)

    return info


def retr_llik(gdat, modlcnts):
    
    if gdat.liketype == 'pois':
        llik = gdat.datacnts * log(modlcnts) - modlcnts
    if gdat.liketype == 'gaus':
        llik = -0.5 * (gdat.datacnts - modlcnts)**2 / gdat.datacnts
    
    return llik


def retr_fluxbias(gdat, spec, indxenerthis):

    # convenience variables
    numbpnts = spec.size
    
    # tolerance factor at the minimum flux
    factlowr = 5. * ones(numbpnts)

    # tolerance factor at the maximum flux
    factuppr = 1.1 * ones(numbpnts)

    slop = (log(factuppr) - log(factlowr)) / (log(gdat.minmflux) - log(gdat.maxmflux))
    offs = log(factuppr) + slop * log(gdat.maxmflux)

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

    angl = sqrt(lgalrttr[None, :, :, None]**2 + bgalrttr[None, :, :, None]**2)
    mapssers = retr_sersprof(spec[:, None, None, None], angl, size)

    return mapssers


def retr_sersprof(spec, angl, size, indx=4):
    
    fact = 3459.5
    sersprof = spec / pi / 40320. / (size / fact)**2 * exp(-(angl / (size / fact))**(1. / indx))
    #sersprof = spec / 2. / pi  / size**2 / indx / sp.special.gamma(2. * indx) * exp(-(angl / size)**(1. / indx))
     
    return sersprof



