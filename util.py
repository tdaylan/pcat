# common imports
from __init__ import *

class interp1d_pick:
    """ class wrapper for piecewise linear function
    """
    def __init__(self, xi, yi, **kwargs):
        self.xi = xi
        self.yi = yi
        self.args = kwargs
        self.f = scipy.interpolate.interp1d(xi, yi, **kwargs)

    def __call__(self, xnew):
        return self.f(xnew)

    def __getstate__(self):
        return self.xi, self.yi, self.args

    def __setstate__(self, state):
        self.f = scipy.interpolate.interp1d(state[0], state[1], **state[2])


def retr_psfnwdth(gdat, psfn, frac):

    wdth = np.zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            psfntemp = psfn[i, :, m]
            indxanglgood = np.argsort(psfntemp)
            intpwdth = max(frac * np.amax(psfntemp), np.amin(psfntemp))
            if intpwdth >= np.amin(psfntemp[indxanglgood]) and intpwdth <= np.amax(psfntemp[indxanglgood]):
                wdthtemp = interp1d_pick(psfntemp[indxanglgood], gdat.binsangl[indxanglgood])(intpwdth)
            else:
                wdthtemp = 0.
            wdth[i, m] = wdthtemp
                        
    return wdth


# lensing-related
def samp_lgalbgalfromtmpl(gdat, probtmpl):
    
    indxpixldraw = choice(gdat.indxpixl, p=probtmpl)
    lgal = gdat.lgalgrid[indxpixldraw] + randn(gdat.sizepixl)
    bgal = gdat.bgalgrid[indxpixldraw] + randn(gdat.sizepixl)
    
    return lgal, bgal


## custom random variables, pdfs, cdfs and icdfs
def cdfn_powr(flux, minm, maxm, slop):
        
    unit = (flux**(1. - slop) - minm**(1. - slop)) / (maxm**(1. - slop) - minm**(1. - slop))
    
    return unit


def cdfn_igam(xdat, slop, cutf):
    
    cdfn = sp.stats.invgamma.cdf(xdat, slop - 1., scale=cutf)
    
    return cdfn


def icdf_dpow(unit, minm, maxm, brek, sloplowr, slopuppr):
    
    if np.isscalar(unit):
        unit = np.array([unit])
    
    faca = 1. / (brek**(sloplowr - slopuppr) * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) / (1. - sloplowr) + (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    facb = faca * brek**(sloplowr - slopuppr) / (1. - sloplowr)

    para = np.empty_like(unit)
    cdfnbrek = facb * (brek**(1. - sloplowr) - minm**(1. - sloplowr))
    indxlowr = np.where(unit <= cdfnbrek)[0]
    indxuppr = np.where(unit > cdfnbrek)[0]
    if indxlowr.size > 0:
        para[indxlowr] = (unit[indxlowr] / facb + minm**(1. - sloplowr))**(1. / (1. - sloplowr))
    if indxuppr.size > 0:
        para[indxuppr] = ((1. - slopuppr) * (unit[indxuppr] - cdfnbrek) / faca + brek**(1. - slopuppr))**(1. / (1. - slopuppr))
    
    return para


def cdfn_dpow(para, minm, maxm, brek, sloplowr, slopuppr):
    
    if np.isscalar(para):
        para = np.array([para])
    
    faca = 1. / (brek**(sloplowr - slopuppr) * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) / (1. - sloplowr) + (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    facb = faca * brek**(sloplowr - slopuppr) / (1. - sloplowr)

    cdfn = np.empty_like(para)
    indxlowr = np.where(para <= brek)[0]
    indxuppr = np.where(para > brek)[0]
    
    if indxlowr.size > 0:
        cdfn[indxlowr] = facb * (para[indxlowr]**(1. - sloplowr) - minm**(1. - sloplowr))
    if indxuppr.size > 0:
        cdfnbrek = facb * (brek**(1. - sloplowr) - minm**(1. - sloplowr))
        cdfn[indxuppr] = cdfnbrek + faca / (1. - slopuppr) * (para[indxuppr]**(1. - slopuppr) - brek**(1. - slopuppr))
    
    return cdfn


def icdf_powr(unit, minm, maxm, slop):
    
    para = (unit * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))
    
    return para


def icdf_igam(xdat, slop, cutf):
    
    icdf = sp.stats.invgamma.ppf(xdat, slop - 1., scale=cutf)
    
    return icdf


def cdfn_expo(para, maxm, scal):

    unit = (1. - np.exp(-para / maxm)) / (1. - np.exp(-maxm / scal))

    return unit


def icdf_expo(unit, maxm, scal):

    para = -scal * np.log(1. - unit * (1. - np.exp(-maxm / scal)))

    return para


def pdfn_expo(xdat, maxm, scal):

    if (xdat > maxm).any():
        pdfn = 0.
    else:
        pdfn = 1. / scal / (1. - np.exp(-maxm / scal)) * np.exp(-xdat / scal)

    return pdfn


def icdf_dexp(cdfn, maxm, scal):
    
    if cdfn < 0.5:
        icdf = -icdf_expo(2. * cdfn, maxm, scal)
    else:
        icdf = icdf_expo(2. * (cdfn - 0.5), maxm, scal)
    
    return icdf


def cdfn_dexp(icdf, maxm, scal):
    
    if icdf < 0.:
        cdfn = cdfn_expo(-icdf, maxm, scal)
    else:
        cdfn = cdfn_expo(icdf, maxm, scal)
    
    return cdfn


def pdfn_dexp(xdat, maxm, scal):
    
    pdfn = 0.5 * pdfn_expo(np.fabs(xdat), maxm, scal)

    return pdfn


def pdfn_dpow(xdat, minm, maxm, brek, sloplowr, slopuppr):
    
    if np.isscalar(xdat):
        xdat = np.array([xdat])
    
    faca = 1. / (brek**(sloplowr - slopuppr) * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) / (1. - sloplowr) + (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    facb = faca * brek**(sloplowr - slopuppr) / (1. - sloplowr)
    
    pdfn = np.empty_like(xdat)
    indxlowr = np.where(xdat <= brek)[0]
    indxuppr = np.where(xdat > brek)[0]
    if indxlowr.size > 0:
        pdfn[indxlowr] = faca * brek**(sloplowr - slopuppr) * xdat[indxlowr]**(-sloplowr)
    if indxuppr.size > 0:
        pdfn[indxuppr] = faca * xdat[indxuppr]**(-slopuppr)
    
    return pdfn


def pdfn_powr(xdat, minm, maxm, slop):
  
    norm = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop))
    
    pdfn = norm * xdat**(-slop)
    
    return pdfn


def pdfn_self(xdat, minm, maxm):
    
    pdfn = 1. / (maxm - minm)
    
    return pdfn


def pdfn_logt(xdat, minm, maxm):
    
    pdfn =  1. / (np.log(maxm) - np.log(minm)) / xdat
    
    return pdfn


def pdfn_igam(xdat, slop, cutf):
    
    pdfn = sp.stats.invgamma.pdf(xdat, slop - 1., scale=cutf)
    
    return pdfn


def cdfn_self(para, minmpara, factpara):
    
    cdfn = (para - minmpara) / factpara
    
    return cdfn


def icdf_self(cdfn, minmpara, factpara):
    
    para = factpara * cdfn + minmpara
    
    return para


def cdfn_lnor(para, meanpara, stdvpara):
   
    cdfn = cdfn_gaus(np.log(para), np.log(meanpara), stdvpara)
    
    return cdfn


def icdf_lnor(cdfn, meanpara, stdvpara):
    
    para = np.exp(icdf_gaus(cdfn, np.log(meanpara), stdvpara))

    return para


def pdfn_lnor(xdat, mean, stdv):
    
    pdfn = pdfn_gaus(np.log(xdat), np.log(mean), stdv)

    return pdfn


def cdfn_gaus(para, meanpara, stdvpara):
   
    cdfn = 0.5  * (1. + sp.special.erf((para - meanpara) / np.sqrt(2) / stdvpara))
    
    return cdfn


def icdf_gaus(cdfn, meanpara, stdvpara):
    
    para = meanpara + stdvpara * np.sqrt(2) * sp.special.erfinv(2. * cdfn - 1.)

    return para


def pdfn_gaus(xdat, mean, stdv):
    
    pdfn = 1. / np.sqrt(2. * pi) / stdv * np.exp(-0.5 * ((xdat - mean) / stdv)**2)

    return pdfn


def cdfn_lgau(para, mean, stdv):
    
    cdfn = cdfn_gaus(np.log(para), np.log(mean), stdv)

    return cdfn


def icdf_lgau(cdfn, mean, stdv):
    
    icdf = np.exp(icdf_gaus(cdfn, np.log(mean), stdv))

    return icdf


def pdfn_lgau(xdat, mean, stdv):
    
    pdfn = pdfn_gaus(np.log(xdat), np.log(mean), stdv)

    return pdfn


def cdfn_eerr(para, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    tranpara = (para - meanpara) / stdvpara
    cdfnnormpara = 0.5 * (sp.special.erf(tranpara / np.sqrt(2.)) + 1.)
    cdfn = (cdfnnormpara - cdfnnormminm) / cdfnnormdiff

    return cdfn


def icdf_eerr(cdfn, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    cdfnnormpara = cdfn * cdfnnormdiff + cdfnnormminm
    tranpara = sp.special.erfinv(2. * cdfnnormpara - 1.) * np.sqrt(2)
    para = tranpara * stdvpara + meanpara
   
    return para


def cdfn_logt(para, minmpara, factpara):

    cdfn = np.log(para / minmpara) / factpara

    return cdfn


def icdf_logt(cdfn, minmpara, factpara):
    
    para = np.exp(cdfn * factpara) * minmpara

    return para


def cdfn_atan(para, minmpara, maxmpara):
    
    cdfn = (np.arctan(para) - np.arctan(minmpara)) / (np.arctan(maxmpara) - np.arctan(minmpara))
    
    return cdfn


def icdf_atan(cdfn, minmpara, maxmpara):

    para = tan((np.arctan(maxmpara) - np.arctan(minmpara)) * cdfn + np.arctan(minmpara))
    
    return para


def pdfn_atan(para, minmpara, maxmpara):

    pdfn = 1. / (para**2 + 1.) / (np.arctan(maxmpara) - np.arctan(minmpara))
    
    return pdfn


def cdfn_fixp(gdat, strgmodl, fixp, thisindxfixp):
    
    scalfixp = getattr(gdat, strgmodl + 'scalfixp')[thisindxfixp]
    
    if scalfixp == 'self' or scalfixp == 'logt' or scalfixp == 'atan':
        
        minmfixp = getattr(gdat, strgmodl + 'minmfixp')[thisindxfixp]
        factfixp = getattr(gdat, strgmodl + 'factfixp')[thisindxfixp]

        if scalfixp == 'self':
            fixpunit = cdfn_self(fixp, minmfixp, factfixp)
        elif scalfixp == 'logt':
            fixpunit = cdfn_logt(fixp, minmfixp, factfixp)

        elif scalfixp == 'atan':
            maxmfixp = getattr(gdat, strgmodl + 'maxmfixp')[thisindxfixp]
            fixpunit = cdfn_atan(fixp, minmfixp, maxmfixp)
    
    elif scalfixp == 'gaus' or scalfixp == 'eerr':
        meanfixp = getattr(gdat, strgmodl + 'meanfixp')[thisindxfixp]
        stdvfixp = getattr(gdat, strgmodl + 'stdvfixp')[thisindxfixp]
        if scalfixp == 'eerr':
            cdfnminmfixp = gdat.cdfnminmfixp[thisindxfixp]
            cdfndifffixp = gdat.cdfndifffixp[thisindxfixp]
            fixpunit = cdfn_eerr(fixp, meanfixp, stdvfixp, cdfnminmfixp, cdfndifffixp)
        else:
            fixpunit = cdfn_gaus(fixp, meanfixp, stdvfixp)

    elif scalfixp == 'pois':
        fixpunit = fixp
    
    if gdat.booldiagmode:
        if fixpunit == 0:
            print 'Warning. CDF is zero.'

    return fixpunit


def icdf_samp(gdat, strgmodl, sampunit, indxsampcomp):
    
    # tobechanged
    samp = np.empty_like(sampunit)
    for scaltype in gdat.listscaltype:
        listindxfixpscal = getattr(gdat, strgmodl + 'listindxfixpscal')[scaltype]
        if len(listindxfixpscal) == 0:
            continue
        samp[listindxfixpscal] = icdf_fixp(gdat, strgmodl, sampunit[listindxfixpscal], scaltype, listindxfixpscal)
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    indxcomp = getattr(gdat, strgmodl + 'indxcomp')
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    if indxsampcomp is not None:
        for l in indxpopl:
            for g in indxcomp[l]:
                indxsamptemp = indxsampcomp[liststrgcomp[l][g]][l]
                # even when numbtrap is 0, indxsampcomp may still not be None, because numbelem might be 0
                if indxsamptemp.size == 0:
                    continue
                samp[indxsamptemp] = icdf_trap(gdat, strgmodl, sampunit[indxsamptemp], samp, listscalcomp[l][g], liststrgcomp[l][g], l)
    
                if gdat.booldiagmode:
                    if not np.isfinite(samp[indxsamptemp]).all():
                        print 'listscalcomp[l][g]'
                        print listscalcomp[l][g]
                        print 'liststrgcomp[l][g]'
                        print liststrgcomp[l][g]
                        print 'lg'
                        print l, g
                        print 'indxsamptemp'
                        print indxsamptemp
                        print 'sampunit[indxsamptemp]'
                        print sampunit[indxsamptemp]
                        raise Exception('')

    return samp

    
def icdf_fixp(gdat, strgmodl, fixpunitscal, scaltype, indxfixpscal):
    
    if scaltype == 'self' or scaltype == 'logt' or scaltype == 'atan':
        minmfixp = getattr(gdat, strgmodl + 'minmfixp')[indxfixpscal]
        factfixp = getattr(gdat, strgmodl + 'factfixp')[indxfixpscal]
    if scaltype == 'self':
        fixp = icdf_self(fixpunitscal, minmfixp, factfixp)
    elif scaltype == 'logt':
        fixp = icdf_logt(fixpunitscal, minmfixp, factfixp)
    elif scaltype == 'atan':
        maxmfixp = getattr(gdat, strgmodl + 'maxmfixp')[indxfixpscal]
        fixp = icdf_atan(fixpunitscal, minmfixp, maxmfixp)
    elif scaltype == 'gaus' or scaltype == 'eerr':
        meanfixp = getattr(gdat, strgmodl + 'meanfixp')[indxfixpscal]
        stdvfixp = getattr(gdat, strgmodl + 'stdvfixp')[indxfixpscal]
        if scaltype == 'eerr':
            cdfnminmfixp = getattr(gdat, strgmodl + 'cdfnminmfixp')[indxfixpscal]
            cdfndifffixp = getattr(gdat, strgmodl + 'cdfndifffixp')[indxfixpscal]
            fixp = icdf_eerr(fixpunitscal, meanfixp, stdvfixp, cdfnminmfixp, cdfndifffixp)
        else:
            fixp = icdf_gaus(fixpunitscal, meanfixp, stdvfixp)
    elif scaltype == 'pois':
        fixp = fixpunitscal

    if gdat.booldiagmode:
        if not np.isfinite(fixp).all():
            namefixp = getattr(gdat, strgmodl + 'namefixp')
            print 'scaltype'
            print scaltype
            print 'fixpunitscal'
            print fixpunitscal
            print 'fixp'
            print fixp
            print 'indxfixpscal'
            print indxfixpscal
            print 'namefixp[indxfixpscal]'
            print namefixp[indxfixpscal]
            print
            raise Exception('')

    return fixp


def retr_lprbpois(data, modl):
    
    lprb = data * np.log(modl) - modl - sp.special.gammaln(data + 1)
    
    return lprb
    
        
def icdf_trap(gdat, strgmodl, cdfn, samp, scalcomp, strgcomp, l):
    
    if scalcomp == 'self' or scalcomp == 'expo' or scalcomp == 'powrslop' or scalcomp == 'dpowslopbrek':
        if scalcomp != 'expo':
            minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
        if scalcomp == 'powrslop':
            maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
            distslop = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]
            
            if gdat.booldiagmode:
                if not np.isfinite(distslop):
                    print 'distslop'
                    print distslop
                    raise Exception('')
                if maxm < minm:
                    print 'strgcomp'
                    print strgcomp
                    print 'minm'
                    print minm
                    print 'maxm'
                    print maxm
                    raise Exception('')
            icdf = icdf_powr(cdfn, minm, maxm, distslop)

        elif scalcomp == 'dpowslopbrek':
            maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
            distbrek = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distbrek')[l]]
            distsloplowr = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distsloplowr')[l]]
            distslopuppr = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslopuppr')[l]]
            icdf = icdf_dpow(cdfn, minm, maxm, distbrek, distsloplowr, distslopuppr)
        elif scalcomp == 'expo':
            maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
            sexp = getattr(gdat, strgmodl + strgcomp + 'distsexppop%d' % l)
            icdf = icdf_expo(cdfn, maxm, sexp)
        else:
            fact = getattr(gdat, strgmodl + 'fact' + strgcomp)
            icdf = icdf_self(cdfn, minm, fact)
    
    if scalcomp == 'logt':
        minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
        fact = getattr(gdat, strgmodl + 'fact' + strgcomp)
        icdf = icdf_logt(cdfn, minm, fact)
    
    if scalcomp.startswith('dexp'):
        if scalcomp == 'dexpscal':
            maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
            scal = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distscal')[l]]
            icdf = icdf_dnp.exp(cdfn, maxm, scal)
    if scalcomp == 'lnormeanstdv':
        distmean = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distmean')[l]]
        diststdv = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'diststdv')[l]]
        icdf = icdf_lnor(cdfn, distmean, diststdv)
    if scalcomp == 'igam':
        distslop = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]
        cutf = getattr(gdat, 'cutf' + strgcomp)
        icdf = icdf_igam(cdfn, distslop, cutf)
    if scalcomp == 'gaus':
        distmean = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distmean')[l]]
        diststdv = samp[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'diststdv')[l]]
        icdf = icdf_gaus(cdfn, distmean, diststdv)
    
    if gdat.booldiagmode:
        if np.isscalar(icdf) and not np.isfinite(icdf) or not np.isscalar(icdf) and not np.isfinite(icdf).any():
            print 'cdfn'
            print cdfn
            print 'l'
            print l
            print 'strgcomp'
            print strgcomp
            print 'scalcomp'
            print scalcomp
            if scalcomp == 'powrslop':
                print 'minm'
                print minm
                print 'maxm'
                print maxm
                print 'distslop'
                print distslop
            print 'icdf'
            print icdf
            raise Exception('')

    return icdf


def cdfn_trap(gdat, gdatmodi, icdf, indxpoplthis):
    
    listscalcomp = gdat.fittlistscalcomp[indxpoplthis]
    cdfn = np.empty(gdat.fittnumbcomp[indxpoplthis])
    for k, strgcomp in enumerate(gdat.fittliststrgcomp[indxpoplthis]):
        
        if listscalcomp[k] == 'self' or listscalcomp[k] == 'dexp' or listscalcomp[k] == 'expo' \
                        or listscalcomp[k] == 'powrslop' or listscalcomp[k] == 'dpowslopbrek':
            minm = getattr(gdat, 'fittminm' + strgcomp)
            if listscalcomp[k] == 'powrslop':
                maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                distslop = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[indxpoplthis]]
                cdfn[k] = cdfn_powr(icdf[k], minm, maxm, distslop)
            elif listscalcomp[k] == 'dpowslopbrek':
                maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                brek = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distbrek')[indxpoplthis]]
                distsloplowr = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distsloplowr')[indxpoplthis]]
                distslopuppr = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslopuppr')[indxpoplthis]]
                cdfn[k] = cdfn_dpow(icdf[k], minm, maxm, brek, distsloplowr, distslopuppr)
            else:
                fact = getattr(gdat, 'fittfact' + strgcomp)
                cdfn[k] = cdfn_self(icdf[k], minm, fact)
        if listscalcomp[k] == 'lnormeanstdv':
            distmean = gdatmodi.samp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[indxpoplthis]]
            diststdv = gdatmodi.samp[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[indxpoplthis]]
            cdfn[k] = cdfn_lnor(icdf[k], distmean, distslop)
        if listscalcomp[k] == 'igam':
            distslop = gdatmodi.samp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[indxpoplthis]]
            cutf = getattr(gdat, 'cutf' + strgcomp)
            cdfn[k] = cdfn_igam(icdf[k], distslop, cutf)
        if listscalcomp[k] == 'gaus':
            distmean = gdatmodi.samp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[indxpoplthis]]
            diststdv = gdatmodi.samp[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[indxpoplthis]]
            cdfn[k] = cdfn_gaus(icdf[k], distmean, diststdv)
    
    return cdfn


### update sampler state
def updt_stat(gdat, gdatmodi):
   
    if gdat.verbtype > 1:
        print 'updt_stat()'
    
    # update the sample and the unit sample vectors
    gdatmodi.thislpritotl = gdatmodi.nextlpritotl
    gdatmodi.thislliktotl = gdatmodi.nextlliktotl
    gdatmodi.thislpostotl = gdatmodi.nextlpostotl
    gdatmodi.thissamp[gdatmodi.indxsampmodi] = np.copy(gdatmodi.nextsamp[gdatmodi.indxsampmodi])
    gdatmodi.thissampunit[gdatmodi.indxsampmodi] = np.copy(gdatmodi.nextsampunit[gdatmodi.indxsampmodi])
    if gdatmodi.thisindxproptype > 0:
        if gdat.verbtype > 1:
            print 'gdatmodi.nextindxelemfull'
            print gdatmodi.nextindxelemfull
        gdatmodi.thisindxelemfull = deepcopy(gdatmodi.nextindxelemfull)
        gdatmodi.thisindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.thisindxelemfull, 'fitt')


def initcompfromstat(gdat, gdatmodi, namerefr):
    
    for l in gdat.fittindxpopl:
        for g, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
            minm = getattr(gdat, 'fittminm' + strgcomp)
            maxm = getattr(gdat, 'fittmaxm' + strgcomp)
            try:
                comp = getattr(gdat, namerefr + strgcomp)[l][0, :]
                if gdat.fittlistscalcomp[l][g] == 'self' or gdat.fittlistscalcomp[l][g] == 'logt':
                    fact = getattr(gdat, 'fittfact' + strgcomp)
                    if gdat.fittlistscalcomp[l][g] == 'self':
                        compunit = cdfn_self(comp, minm, fact)
                    if gdat.fittlistscalcomp[l][g] == 'logt':
                        compunit = cdfn_logt(comp, minm, fact)
                if gdat.fittlistscalcomp[l][g] == 'expo':
                    scal = getattr(gdat, 'fittgangdistsexp')
                    maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                    compunit = cdfn_expo(icdf, maxm, scal)
                if gdat.fittlistscalcomp[l][g] == 'powrslop' or gdat.fittlistscalcomp[l][g] == 'igam':
                    slop = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[l]]
                    if gdat.fittlistscalcomp[l][g] == 'powrslop':
                        compunit = cdfn_powr(comp, minm, maxm, slop)
                    if gdat.fittlistscalcomp[l][g] == 'igam':
                        cutf = getattr(gdat, 'cutf' + strgcomp)
                        compunit = cdfn_igam(comp, slop, cutf)
                if gdat.fittlistscalcomp[l][g] == 'dpowslopbrek':
                    brek = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distbrek')[l]]
                    sloplowr = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distsloplowr')[l]]
                    slopuppr = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslopuppr')[l]]
                    compunit = cdfn_powr(comp, minm, maxm, brek, sloplowr, slopuppr)
                if gdat.fittlistscalcomp[l][g] == 'gaus':
                    distmean = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[l]]
                    diststdv = gdatmodi.thissamp[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[l]]
                    compunit = cdfn_gaus(comp, distmean, diststdv)
            except:
                if gdat.verbtype > 0:
                    print 'Initialization from the reference catalog failed for %s. Sampling randomly...' % strgcomp
                compunit = rand(gdatmodi.thissamp[gdat.fittindxfixpnumbelem[l]].astype(int))
            gdatmodi.thissampunit[gdatmodi.thisindxsampcomp[strgcomp][l]] = compunit


### find the spectra of sources
def retr_spec(gdat, flux, sind=None, curv=None, expc=None, sindcolr=None, elin=None, edisintp=None, sigm=None, gamm=None, spectype='powr', plot=False):
    
    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if plot:
            meanener = gdat.meanenerplot
        else:
            meanener = gdat.meanener

        if spectype == 'gaus':
            spec = 1. / edis[None, :] / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanener[:, None] - elin[None, :]) / edis[None, :])**2)
        if spectype == 'voig':
            args = (gdat.meanener[:, None] + 1j * gamm[None, :]) / np.sqrt(2.) / sigm[None, :]
            spec = 1. / sigm[None, :] / np.sqrt(2. * pi) * flux[None, :] * real(scipy.special.wofz(args))
        if spectype == 'edis':
            edis = edisintp(elin)[None, :]
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanener[:, None] - elin[None, :]) / edis)**2)
        if spectype == 'pvoi':
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanener[:, None] - elin[None, :]) / edis)**2)
        if spectype == 'lore':
            spec = 1. / edis / np.sqrt(2. * pi) * flux[None, :] * np.exp(-0.5 * ((gdat.meanener[:, None] - elin[None, :]) / edis)**2)
        if spectype == 'powr':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :])
        if spectype == 'colr':
            if plot:
                spec = np.zeros((gdat.numbenerplot, flux.size))
            else:
                spec = np.empty((gdat.numbener, flux.size))
                for i in gdat.indxener:
                    if i < gdat.indxenerpivt:
                        spec[i, :] = flux * (gdat.meanener[i] / gdat.enerpivt)**(-sindcolr[i])
                    elif i == gdat.indxenerpivt:
                        spec[i, :] =  flux
                    else:
                        spec[i, :] = flux * (gdat.meanener[i] / gdat.enerpivt)**(-sindcolr[i-1])
        if spectype == 'curv':
            spec = flux[None, :] * meanener[:, None]**(-sind[None, :] - gdat.factlogtenerpivt[:, None] * curv[None, :])
        if spectype == 'expc':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :]) * np.exp(-(meanener - gdat.enerpivt)[:, None] / expc[None, :])
    
    return spec


### find the surface brightness due to one point source
def retr_sbrtpnts(gdat, lgal, bgal, spec, psfnintp, indxpixlelem):
    
    # calculate the distance to all pixels from each point source
    dist = retr_angldistunit(gdat, lgal, bgal, indxpixlelem)
    
    # interpolate the PSF onto the pixels
    if gdat.kernevaltype == 'ulip':
        psfntemp = psfnintp(dist)
    if gdat.kernevaltype == 'bspx':
        pass

    # scale by the PS spectrum
    sbrtpnts = spec[:, None, None] * psfntemp
    
    return sbrtpnts


### find the set of pixels in proximity to a position on the map
def retr_indxpixlelemconc(gdat, strgmodl, dictelem, l):

    elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    boolelemlght = getattr(gdat, strgmodl + 'boolelemlght')
    
    lgal = dictelem[l]['lgal']
    bgal = dictelem[l]['bgal']
    
    if boolelemlght[l]:
        varbampl = abs(dictelem[l]['spec'][gdat.indxenerpivt, :])
    if elemtype[l] == 'lens':
        varbampl = dictelem[l]['defs']
    if elemtype[l].startswith('clus'):
        varbampl = dictelem[l]['nobj']
    
    if elemspatevaltype[l] == 'locl':
        listindxpixlelem = [[] for k in range(lgal.size)]
        for k in range(lgal.size):
            indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
            
            indxfluxproxtemp = np.digitize(varbampl[k], gdat.binsprox)
            if indxfluxproxtemp > 0:
                indxfluxproxtemp -= 1
            if indxfluxproxtemp == gdat.binsprox.size - 1:
                print 'Warning! Index of the proximity pixel list overflew. Taking the largest list...'
                print 'varbampl[k]'
                print varbampl[k]
                print 'gdat.binsprox'
                print gdat.binsprox
                print
                indxfluxproxtemp -= 1
            if gdat.verbtype > 1:
                print 'k'
                print k
                print 'indxpixlpnts'
                print indxpixlpnts
                print 'indxfluxproxtemp'
                print indxfluxproxtemp
            
            indxpixlelem = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
            
            if isinstance(indxpixlelem, int):
                indxpixlelem = gdat.indxpixl
            
            listindxpixlelem[k] = indxpixlelem
            if gdat.verbtype > 1:
                print 'dictelem[l][lgal][k]'
                print dictelem[l]['lgal'][k]
                print 'dictelem[l][bgal][k]'
                print dictelem[l]['bgal'][k]
                print 'varbampl'
                print varbampl
                print 'indxpixlelem'
                summgene(indxpixlelem)
                print

        listindxpixlelemconc = np.unique(np.concatenate(listindxpixlelem))
    else:
        listindxpixlelemconc = gdat.indxpixl
        listindxpixlelem = gdat.indxpixl
    
    return listindxpixlelem, listindxpixlelemconc


### find the distance between two points on the map
def retr_angldistunit(gdat, lgal, bgal, indxpixlelem, retranglcosi=False):
   
    if gdat.pixltype == 'heal':
        xdat, ydat, zaxi = retr_unit(lgal, bgal)
        anglcosi = gdat.xdatgrid[indxpixlelem] * xdat + gdat.ydatgrid[indxpixlelem] * ydat + gdat.zaxigrid[indxpixlelem] * zaxi
        
        if retranglcosi:
            return anglcosi
        else:
            angldist = np.arccos(anglcosi)
            return angldist
    
    else:
        angldist = np.sqrt((lgal - gdat.lgalgrid[indxpixlelem])**2 + (bgal - gdat.bgalgrid[indxpixlelem])**2)
        
        return angldist
    

### find the pixel index of a point on the map
def retr_indxpixl(gdat, bgal, lgal):

    if gdat.pixltype == 'heal':
        indxpixl = gdat.pixlcnvt[ang2pix(gdat.numbsideheal, np.pi / 2. - bgal, lgal)]
        if gdat.booldiagmode:
            if (indxpixl == -1).any():  
                raise Exception('pixlcnvt went negative!')

    if gdat.pixltype == 'cart':
        indxlgcr = np.floor(gdat.numbsidecart * (lgal - gdat.minmlgaldata) / 2. / gdat.maxmgangdata).astype(int)
        indxbgcr = np.floor(gdat.numbsidecart * (bgal - gdat.minmbgaldata) / 2. / gdat.maxmgangdata).astype(int)

        if np.isscalar(indxlgcr):
            if indxlgcr < 0:
                indxlgcr = 0
            if indxlgcr >= gdat.numbsidecart:
                indxlgcr = gdat.numbsidecart - 1
        else:
            indxlgcr[np.where(indxlgcr < 0)] = 0
            indxlgcr[np.where(indxlgcr >= gdat.numbsidecart)] = gdat.numbsidecart - 1
            
        if np.isscalar(indxbgcr):
            if indxbgcr < 0:
                indxbgcr = 0
            if indxbgcr >= gdat.numbsidecart:
                indxbgcr = gdat.numbsidecart - 1
        else:
            indxbgcr[np.where(indxbgcr < 0)] = 0
            indxbgcr[np.where(indxbgcr >= gdat.numbsidecart)] = gdat.numbsidecart - 1
            
        indxpixl = indxlgcr * gdat.numbsidecart + indxbgcr
    
    # convert to an index of non-zero exposure pixels
    #indxpixl = gdat.indxpixlroficnvt[indxpixl]

    return indxpixl


## obtain count maps
def retr_cntp(gdat, sbrt):
   
    cntp = sbrt * gdat.expo * gdat.apix
    if gdat.enerdiff:
        cntp *= gdat.deltener[:, None, None] 
        
    return cntp


## plotting
### construct path for plots
def retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, strgplot, nameinte=''):
    
    if strgmodl == 'true' or strgstat == '':
        path = gdat.pathinit + nameinte + strgplot + '.pdf'
    elif strgstat == 'pdfn' or strgstat == 'mlik':
        path = gdat.pathplotrtag + strgpdfn + '/finl/' + nameinte + strgstat + strgplot + '.pdf'
    elif strgstat == 'this':
        path = gdat.pathplotrtag + strgpdfn + '/fram/' + nameinte + strgstat + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


### determine the marker size
def retr_mrkrsize(gdat, compampl, namefeatampl):
    
    minm = getattr(gdat, 'minm' + namefeatampl) 
    maxm = getattr(gdat, 'maxm' + namefeatampl)
    mrkrsize = (np.sqrt(compampl) - np.sqrt(minm)) / (np.sqrt(maxm) - np.sqrt(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


## experiment specific
def retr_psfphubb(gdat):

    # temp
    gdat.psfpexpr = np.array([0.080, 0.087]) / gdat.anglfact


def retr_psfpchan(gdat):

    # temp
    #gdat.psfpexpr = np.array([0.25, 0.3, 0.4, 0.6, 0.7]) / gdat.anglfact
    if gdat.numbenerfull == 5:
        gdat.psfpexpr = np.array([0.424 / gdat.anglfact, 2.75, 0.424 / gdat.anglfact, 2.59, 0.440 / gdat.anglfact, 2.47, 0.457 / gdat.anglfact, 2.45, 0.529 / gdat.anglfact, 3.72])
    if gdat.numbenerfull == 2:
        gdat.psfpexpr = np.array([0.427 / gdat.anglfact, 2.57, 0.449 / gdat.anglfact, 2.49])
    #gdat.psfpchan = gdat.psfpexpr[(2 * gdat.indxenerincl[:, None] + np.arange(2)[None, :]).flatten()] 
    #gdat.psfpexpr = np.array([0.25 / gdat.anglfact, 
    #                       0.30 / gdat.anglfacti\
    #                       0.40 / gdat.anglfacti\
    #                       0.60 / gdat.anglfacti\
    #                       0.70 / gdat.anglfacti
    #gdat.psfpexpr = np.array([0.35 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.e-1, 2.])
    #gdat.psfpexpr = np.array([0.25 / gdat.anglfact, 2.0e-1, 1.9, \
    #                       0.30 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.40 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.60 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.70 / gdat.anglfact, 1.0e-1, 2.0])
   

def retr_psfpsdyn(gdat):

    gdat.psfpexpr = np.array([0.05 / gdat.anglfact])
   

def retr_psfpferm(gdat):
   
    if gdat.anlytype.startswith('rec8'):
        path = gdat.pathdata + 'expr/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    else:
        path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
    irfn = astropy.io.fits.getdata(path, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = np.sqrt(minmener * maxmener)

    numbpsfpscal = 3
    numbpsfpform = 5
    
    fermscal = np.zeros((gdat.numbevtt, numbpsfpscal))
    fermform = np.zeros((gdat.numbener, gdat.numbevtt, numbpsfpform))
    
    strgpara = ['score', 'gcore', 'stail', 'gtail', 'ntail']
    for m in gdat.indxevtt:
        if gdat.anlytype.startswith('rec8'):
            irfn = astropy.io.fits.getdata(path, 1 + 3 * gdat.indxevttincl[m])
            fermscal[m, :] = astropy.io.fits.getdata(path, 2 + 3 * gdat.indxevttincl[m])['PSFSCALE']
        else:
            if m == 1:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_front.fits'
            elif m == 0:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
            else:
                continue
            irfn = astropy.io.fits.getdata(path, 1)
            fermscal[m, :] = astropy.io.fits.getdata(path, 2)['PSFSCALE']
        for k in range(numbpsfpform):
            fermform[:, m, k] = interp1d_pick(enerirfn, np.mean(irfn[strgpara[k]].squeeze(), axis=0))(gdat.meanener)
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 4] = 1. / (1. + fermform[i, m, 4] * fermform[i, m, 2]**2 / fermform[i, m, 0]**2)

    # calculate the scale factor
    gdat.fermscalfact = np.sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # store the fermi PSF parameters
    gdat.psfpexpr = np.zeros(gdat.numbener * gdat.numbevtt * numbpsfpform)
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k]
    

def retr_refrchaninit(gdat):
    
    gdat.indxrefr = np.arange(gdat.numbrefr)
    
    gdat.dictrefr = []
    for q in gdat.indxrefr:
        gdat.dictrefr.append(dict())
    
    gdat.namefeatsignrefr = 'flux'
    
    gdat.refrlegdelem = ['Xue+2011', 'Wolf+2008']
    
    gdat.refrlistnamefeatampl[0] = 'flux'
    gdat.refrlistnamefeatampl[1] = 'magt'
    gdat.listnamerefr += ['xu11', 'wo08']
    
    setattr(gdat, 'minmotyp', 0.)
    setattr(gdat, 'maxmotyp', 1.)
    setattr(gdat, 'lablotyp', 'O')
    setattr(gdat, 'factotypplot', 1.)
    setattr(gdat, 'scalotypplot', 'self')
    
    setattr(gdat, 'lablotypxu11', 'O')
    for name in gdat.listnamerefr:
        setattr(gdat, 'minmotyp' + name, 0.)
        setattr(gdat, 'maxmotyp' + name, 1.)
    
    if False and gdat.strgcnfg == 'pcat_chan_inpt_home4msc':
        print 'ECDFS_Cross_ID_Hsu2014'
        with open(gdat.pathinpt + 'ECDFS_Cross_ID_Hsu2014.txt', 'r') as thisfile:
            for k, line in enumerate(thisfile):
                if k < 18:
                    continue
                rasccand =line[2]
                declcand =line[2]
       
    gdat.refrliststrgfeat[0] += ['lgal', 'bgal', 'flux', 'sind', 'otyp', 'lumi']
    gdat.refrliststrgfeat[1] += ['lgal', 'bgal', 'magt', 'reds', 'otyp']


def retr_refrchanfinl(gdat):
    
    booltemp = False
    if gdat.anlytype.startswith('extr'):
        if gdat.numbsidecart == 300:
            gdat.numbpixllgalshft[0] = 1490
            gdat.numbpixlbgalshft[0] = 1430
        else:
            booltemp = True
    elif gdat.anlytype.startswith('home'):
        gdat.numbpixllgalshft[0] = 0
        gdat.numbpixlbgalshft[0] = 0
    
        if gdat.numbsidecart == 600:
            pass
        elif gdat.numbsidecart == 100:
            indxtile = int(gdat.anlytype[-4:])
            numbsidecntr = int(gdat.anlytype[8:12])
            numbtileside = numbsidecntr / gdat.numbsidecart
            indxtilexaxi = indxtile // numbtileside
            indxtileyaxi = indxtile % numbtileside
            gdat.numbpixllgalshft[0] += indxtilexaxi * gdat.numbsidecart
            gdat.numbpixlbgalshft[0] += indxtileyaxi * gdat.numbsidecart
        elif gdat.numbsidecart == 300:
            gdat.numbpixllgalshft[0] += 150
            gdat.numbpixlbgalshft[0] += 150
        else:
            booltemp = True
    else:
        booltemp = True

    if booltemp:
        print 'gdat.numbsidecart'
        print gdat.numbsidecart
        raise Exception('Reference elements cannot be aligned with the spatial axes!')
    
    ## WCS object for rotating reference elements into the ROI
    if gdat.numbener == 2:
        gdat.listpathwcss[0] = gdat.pathinpt + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
    else:
        gdat.listpathwcss[0] = gdat.pathinpt + '0.5-0.91028_flux_%sMs.img' % gdat.anlytype[4]
    
    # Xue et al. (2011)
    #with open(gdat.pathinpt + 'chancatl.txt', 'r') as thisfile:
    pathfile = gdat.pathinpt + 'Xue2011.fits'
    hdun = pf.open(pathfile)
    hdun.info()
    #print repr(hdun[1].header)
    lgalchan = hdun[1].data['_Glon'] / 180. * pi
    bgalchan = hdun[1].data['_Glat'] / 180. * pi
    fluxchansoft = hdun[1].data['SFlux']
    fluxchanhard = hdun[1].data['HFlux']
    objttypechan = hdun[1].data['Otype']
    gdat.refrlumi[0][0] = hdun[1].data['Lx']
    
    # position
    gdat.refrlgal[0] = lgalchan
    gdat.refrbgal[0] = bgalchan

    # spectra
    gdat.refrspec = [[np.zeros((3, gdat.numbener, lgalchan.size))]]
    if gdat.numbener == 2:
        gdat.refrspec[0][0, 0, :] = fluxchansoft * 0.624e9
        gdat.refrspec[0][0, 1, :] = fluxchanhard * 0.624e9 / 16.
    else:
        gdat.refrspec[0][0, :, :] = 2. * fluxchansoft[None, :] * 0.624e9
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :]
   
    # fluxes
    gdat.refrflux[0] = gdat.refrspec[0][:, gdat.indxenerpivt, :]

    # spectral indices
    if gdat.numbener > 1:
        gdat.refrsind[0] = -np.log(gdat.refrspec[0][0, 1, :] / gdat.refrspec[0][0, 0, :]) / np.log(np.sqrt(7. / 2.) / np.sqrt(0.5 * 2.))

    ## object type
    objttypechantemp = np.zeros(lgalchan.size) - 1.
    indx = np.where(objttypechan == 'AGN')[0]
    objttypechantemp[indx] = 0.165
    indx = np.where(objttypechan == 'Galaxy')[0]
    objttypechantemp[indx] = 0.495
    indx = np.where(objttypechan == 'Star')[0]
    objttypechantemp[indx] = 0.835
    gdat.refrotyp[0][0] = objttypechantemp

    # Wolf et al. (2011)
    path = gdat.pathdata + 'inpt/Wolf2008.fits'
    data = astropy.io.fits.getdata(path)
    gdat.refrlgal[1] = np.deg2rad(data['_Glon'])
    gdat.refrlgal[1] = ((gdat.refrlgal[1] - pi) % (2. * pi)) - pi
    gdat.refrbgal[1] = np.deg2rad(data['_Glat'])
    gdat.refrmagt[1][0] = data['Rmag']
    gdat.refrreds[1][0] = data['MCz']
  
    #listname = []
    #for k in range(data['MCclass'].size):
    #    if not data['MCclass'][k] in listname:
    #        listname.append(data['MCclass'][k])
    listname = ['Galaxy', 'Galaxy  (Uncl!)', 'QSO     (Gal?)', 'Galaxy  (Star?)', 'Star', 'Strange Object', 'QSO', 'WDwarf']
    gdat.refrotyp[1][0] = np.zeros_like(gdat.refrreds[1][0]) - 1. 
    for k, name in enumerate(listname):
        indx = np.where(data['MCclass'] == name)[0]
        gdat.refrotyp[1][0][indx] = k / 10.
    
    # error budget
    for name in ['lgal', 'bgal', 'sind', 'otyp', 'lumi', 'magt', 'reds']:
        refrtile = [[] for q in gdat.indxrefr]
        refrfeat = getattr(gdat, 'refr' + name)
        for q in gdat.indxrefr:
            if len(refrfeat[q]) > 0:
                refrtile[q] = np.tile(refrfeat[q], (3, 1))
        setattr(gdat, 'refr' + name, refrtile)
        

def retr_refrferminit(gdat):
    
    gdat.listnamerefr += ['ac15', 'ma05']
    gdat.indxrefr = np.arange(gdat.numbrefr)
    
    gdat.refrlegdelem = ['Acero+2015', 'Manchester+2005']

    gdat.refrlistnamefeatampl[0] = 'flux'
    gdat.refrlistnamefeatampl[1] = 'flux0400'
    gdat.namefeatsignrefr = 'flux'
    
    setattr(gdat, 'lablcurvac15', '%s_{3FGL}' % gdat.lablcurv)
    setattr(gdat, 'lablexpcac15', 'E_{c,3FGL}')
    
    for name in gdat.listnamerefr:
        setattr(gdat, 'minmcurv' + name, -1.)
        setattr(gdat, 'maxmcurv' + name, 1.)
        setattr(gdat, 'minmexpc' + name, 0.1)
        setattr(gdat, 'maxmexpc' + name, 10.)
   
    gdat.refrliststrgfeat[0] += ['lgal', 'bgal', 'flux', 'sind', 'curv', 'expc', 'tvar', 'etag', 'styp', 'sindcolr0001', 'sindcolr0002']
    gdat.refrliststrgfeat[1] += ['lgal', 'bgal', 'flux0400', 'per0', 'per1']


def retr_refrfermfinl(gdat):

    gdat.minmstyp = -0.5
    gdat.maxmstyp = 3.5
    gdat.lablstyp = 'S'
    gdat.factstypplot = 1.
    gdat.scalstypplot = 'self'
    
    gdat.minmtvar = 0.
    gdat.maxmtvar = 400.
    gdat.labltvar = 'T'
    gdat.facttvarplot = 1.
    gdat.scaltvarplot = 'logt'
    
    # Acero+2015
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'
    fgl3 = astropy.io.fits.getdata(path)
    
    gdat.refrlgal[0] = np.deg2rad(fgl3['glon'])
    gdat.refrlgal[0] = np.pi - ((gdat.refrlgal[0] - np.pi) % (2. * np.pi))
    gdat.refrbgal[0] = np.deg2rad(fgl3['glat'])
    
    gdat.refrnumbelemfull = gdat.refrlgal[0].size

    gdat.refrspec = [np.empty((3, gdat.numbener, gdat.refrlgal[0].size))]
    gdat.refrspec[0][0, :, :] = np.stack((fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000']))[gdat.indxenerincl, :] / gdat.deltener[:, None]
    
    fgl3specstdvtemp = np.stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.refrspec[0][np.where(np.isfinite(gdat.refrspec[0]) == False)] = 0.
    
    gdat.refrflux[0] = gdat.refrspec[0][:, gdat.indxenerpivt, :]
    gdat.refrsindcolr0001[0] = -np.log(gdat.refrspec[0][:, 1, :] / gdat.refrflux[0]) / np.log(gdat.meanener[1] / gdat.enerpivt)
    gdat.refrsindcolr0002[0] = -np.log(gdat.refrspec[0][:, 2, :] / gdat.refrflux[0]) / np.log(gdat.meanener[2] / gdat.enerpivt)
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = np.deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(np.cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(np.sin(fgl3anglstdv))

    gdat.refretag[0] = np.zeros(gdat.refrlgal[0].size, dtype=object)
    for k in range(gdat.refrlgal[0].size):
        gdat.refretag[0][k] = '%s, %s, %s' % (fgl3['Source_Name'][k], fgl3['CLASS1'][k], fgl3['ASSOC1'][k])
    gdat.refrtvar[0] = fgl3['Variability_Index']
    
    gdat.refrstyp[0] = np.zeros_like(gdat.refrlgal[0]) - 1
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PowerLaw        ')] = 0
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'LogParabola     ')] = 1
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PLExpCutoff     ')] = 2
    gdat.refrstyp[0][np.where(fgl3['SpectrumType'] == 'PLSuperExpCutoff')] = 3
    indx = np.where(gdat.refrstyp[0] == -1)[0]
    if indx.size > 0:
        print 'indx'
        print indx
        print 'fgl3[SpectrumType][indx]'
        print fgl3['SpectrumType'][indx]
        raise Exception('')
    gdat.refrsind[0] = fgl3['Spectral_Index']
    gdat.refrcurv[0] = fgl3['beta']
    gdat.refrexpc[0] = fgl3['Cutoff'] * 1e-3
    
    gdat.refrcurv[0][np.where(np.logical_not(np.isfinite(gdat.refrcurv[0])))] = -10.
    gdat.refrexpc[0][np.where(np.logical_not(np.isfinite(gdat.refrexpc[0])))] = 0.
    
    gdat.refrsind[0] = np.tile(gdat.refrsind[0], (3, 1)) 
    gdat.refrcurv[0] = np.tile(gdat.refrcurv[0], (3, 1)) 
    gdat.refrexpc[0] = np.tile(gdat.refrexpc[0], (3, 1)) 

    # Manchester+2005
    path = gdat.pathdata + 'inpt/Manchester2005.fits'
    data = astropy.io.fits.getdata(path)
   
    gdat.refrlgal[1] = np.deg2rad(data['glon'])
    gdat.refrlgal[1] = ((gdat.refrlgal[1] - np.pi) % (2. * np.pi)) - np.pi
    gdat.refrbgal[1] = np.deg2rad(data['glat'])
    
    gdat.refrper0[1] = data['P0']
    gdat.refrper1[1] = data['P1']
    gdat.refrflux0400[1] = data['S400']
    #gdat.refrdism[1] = data['DM']
    #gdat.refrdlos[1] = data['Dist']

    # error budget
    for name in ['lgal', 'bgal', 'per0', 'per1', 'flux0400', 'tvar', 'styp']:
        refrtile = [[] for q in gdat.indxrefr]
        refrfeat = getattr(gdat, 'refr' + name)
        for q in gdat.indxrefr:
            if len(refrfeat[q]) > 0:
                refrtile[q] = np.tile(refrfeat[q], (3, 1))
        setattr(gdat, 'refr' + name, refrtile)


def retr_singgaus(scaldevi, sigc):
    
    psfn = 1. / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_singking(scaldevi, sigc, gamc):
    
    psfn = 1. / 2. / np.pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc)

    return psfn


def retr_doubgaus(scaldevi, frac, sigc, sigt):
    
    psfn = frac / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2) + (1. - frac) / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_gausking(scaldevi, frac, sigc, sigt, gamt):

    psfn = frac / 2. / np.pi / sigc**2 * np.exp(-0.5 * scaldevi**2 / sigc**2) + (1. - frac) / 2. / np.pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_doubking(scaldevi, frac, sigc, gamc, sigt, gamt):

    psfn = frac / 2. / np.pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc) + \
                            (1. - frac) / 2. / np.pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_lgalbgal(gang, aang):
    
    lgal = gang * np.cos(aang)
    bgal = gang * np.sin(aang)

    return lgal, bgal


def retr_gang(lgal, bgal):
    
    gang = np.arccos(np.cos(lgal) * np.cos(bgal))

    return gang


def retr_aang(lgal, bgal):

    aang = np.arctan2(bgal, lgal)

    return aang


def show_samp(gdat, gdatmodi, strgstat='this', strgmodl='fitt'):
    
    namepara = getattr(gdat, strgmodl + 'namepara')
    scalpara = getattr(gdat, strgmodl + 'scalpara')
    indxpara = getattr(gdat, strgmodl + 'indxpara')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    
    strgpfix = retr_strgpfix(strgstat, strgmodl)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl)
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        indxsampcomp = getattr(gdatobjt, strgpfix + 'indxsampcomp')
    sampunit = getattr(gdatobjt, strgpfix + 'sampunit')
    samp = getattr(gdatobjt, strgpfix + 'samp')

    print 'strgstat: ' + strgstat
    print 'strgmodl: ' + strgmodl
    print '%5s %20s %30s %30s %15s' % ('index', 'namepara', 'sampunit', 'samp', 'scalpara')
    for k in indxpara:
        if numbtrap > 0:
            if k in np.concatenate(indxsampcomp['lgal']):
                print
        print '%5d %20s %30g %30g %15s' % (k, namepara[k], sampunit[k], samp[k], scalpara[k])
    

def prop_stat(gdat, gdatmodi, strgmodl, thisindxelem=None, thisindxpopl=None, brth=False, deth=False):
 
    if gdat.verbtype > 1:
        print 'prop_stat()'
    
    #indxproptype
    # within, birth, death, split, merge
    # 0, 1, 2, 3, 4

    elemtype = getattr(gdat, strgmodl + 'elemtype')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        numbcomp = getattr(gdat, strgmodl + 'numbcomp')
        indxcomp = getattr(gdat, strgmodl + 'indxcomp')
        minmnumbelem = getattr(gdat, strgmodl + 'minmnumbelem')
        maxmnumbelem = getattr(gdat, strgmodl + 'maxmnumbelem')
    namepara = getattr(gdat, strgmodl + 'namepara')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')
    if psfnevaltype != 'none':
        indxfixppsfp = getattr(gdat, strgmodl + 'indxfixppsfp')
    indxfixpbacp = getattr(gdat, strgmodl + 'indxfixpbacp')
    
    indxfixphostlght = getattr(gdat, strgmodl + 'indxfixphostlght')
    indxfixphostlens = getattr(gdat, strgmodl + 'indxfixphostlens')
    indxfixpsour = getattr(gdat, strgmodl + 'indxfixpsour')
    indxfixpextr = getattr(gdat, strgmodl + 'indxfixpextr')
    indxfixplenp = getattr(gdat, strgmodl + 'indxfixplenp')
    
    if numbtrap > 0:
        indxfixpdist = getattr(gdat, strgmodl + 'indxfixpdist')
        liststrgfeatdefa = getattr(gdat, strgmodl + 'liststrgfeatdefa')
        indxfixpmeanelem = getattr(gdat, strgmodl + 'indxfixpmeanelem')
        indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
        maxmnumbelempopl = getattr(gdat, strgmodl + 'maxmnumbelempopl')
        spatdisttype = getattr(gdat, strgmodl + 'spatdisttype')
        boolelempsfn = getattr(gdat, strgmodl + 'boolelempsfn')
        boolelemsbrtdfnc = getattr(gdat, strgmodl + 'boolelemsbrtdfnc')
        boolcompposi = getattr(gdat, strgmodl + 'boolcompposi')
        indxcompampl = getattr(gdat, strgmodl + 'indxcompampl')
    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
    indxsampcompinit = getattr(gdat, strgmodl + 'indxsampcompinit')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    convdiffanyy = getattr(gdat, strgmodl + 'convdiffanyy')
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    numbtrappoplcuml = getattr(gdat, strgmodl + 'numbtrappoplcuml')
    
    spectype = getattr(gdat, strgmodl + 'spectype')
        
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, 'this', strgmodl)

    strgpfixthis = retr_strgpfix('this', strgmodl)
    strgpfixnext = retr_strgpfix('next', strgmodl)

    thissamp = getattr(gdatobjt, strgpfixthis + 'samp')
    thissampunit = getattr(gdatobjt, strgpfixthis + 'sampunit')
    if numbtrap > 0:
        boolelemlght = getattr(gdat, strgmodl + 'boolelemlght')
        thisindxelemfull = getattr(gdatobjt, strgpfixthis + 'indxelemfull')
        if gdat.booldiagmode:
            for l in gdat.fittindxpopl:
                if len(thisindxelemfull[l]) > len(set(thisindxelemfull[l])):
                    raise Exception('Repeating entry in the element index list!')

        thisindxsampcomp = retr_indxsampcomp(gdat, thisindxelemfull, strgmodl)
        setattr(gdatobjt, strgpfixthis + 'indxsampcomp', thisindxsampcomp)
    else:
        thisindxsampcomp = None
    
    gdatmodi.thisboolpropfilt = True 

    # index of the population in which a transdimensional proposal will be attempted
    if numbtrap > 0:
        if thisindxpopl is None:
            gdatmodi.indxpopltran = choice(indxpopl)
        else:
            gdatmodi.indxpopltran = thisindxpopl
        numbelemtemp = thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]]
        maxmnumbelemtemp = maxmnumbelem[gdatmodi.indxpopltran]
        minmnumbelemtemp = minmnumbelem[gdatmodi.indxpopltran]
    
    # forced death or birth does not check for the prior on the dimensionality on purpose!
    if numbtrap > 0 and (deth or brth or rand() < gdat.probtran) and not (numbelemtemp == minmnumbelemtemp and numbelemtemp == maxmnumbelemtemp):

        if brth or deth or rand() < gdat.probbrde or numbelemtemp == maxmnumbelemtemp and numbelemtemp == 1 or numbelemtemp == 0:
            
            ## births and deaths
            if numbelemtemp == maxmnumbelemtemp or deth:
                gdatmodi.thisindxproptype = 2
            elif numbelemtemp == minmnumbelemtemp or brth:
                gdatmodi.thisindxproptype = 1
            else:
                if rand() < 0.5:
                    gdatmodi.thisindxproptype = 1
                else:
                    gdatmodi.thisindxproptype = 2

        else:
            ## splits and merges
            if numbelemtemp == minmnumbelemtemp or numbelemtemp < 2:
                gdatmodi.thisindxproptype = 3
            elif numbelemtemp == maxmnumbelemtemp:
                gdatmodi.thisindxproptype = 4
            else:
                if rand() < 0.5:
                    gdatmodi.thisindxproptype = 3
                else:
                    gdatmodi.thisindxproptype = 4
    else:
        
        if gdatmodi.boolburn and gdat.boolevolstdp:
            gdatmodi.stdpprev = np.copy(gdatmodi.stdp)
            
            # increase all proposal sizes
            gdatmodi.stdp += gdatmodi.stdp * 0.1 * np.random.random(gdat.numbstdp)
            # decrease those proposal sizes that inreased in the last rejected step
            if gdatmodi.indxstdpincr.size > 0:
                gdatmodi.stdp[gdatmodi.indxstdpincr] -= gdatmodi.stdp[gdatmodi.indxstdpincr] * 0.1 * np.random.randn(gdatmodi.indxstdpincr.size)
        else:
            stdp = gdatmodi.stdp
        
        if gdat.booldiagmode and (gdatmodi.stdp > 1e2).any():
            print 'Warning! gdatmodi.stdp went too large.'
            print 'gdatmodi.stdp'
            print gdatmodi.stdp
            summgene(gdatmodi.stdp)
            #raise Exception('')

        # get the indices of the current parameter vector
        if numbtrap > 0:
            thisindxsampfull = np.concatenate([gdat.fittindxfixpstdv] + thisindxsampcomp['comp'])
        else:
            thisindxsampfull = gdat.fittindxfixpstdv
        
        thisstdp = gdatmodi.stdp[gdat.indxstdppara[thisindxsampfull]]
        if not np.isfinite(thisstdp).all():
            print 'thisindxsampfull'
            summgene(thisindxsampfull)
            print 'gdat.indxstdppara'
            summgene(gdat.indxstdppara)
            print 'gdat.indxstdppara[thisindxsampfull]'
            summgene(gdat.indxstdppara[thisindxsampfull])
            print 'gdatmodi.stdp'
            summgene(gdatmodi.stdp)
            print 'gdatmodi.stdp[gdat.indxstdppara[thisindxsampfull]]'
            print gdatmodi.stdp[gdat.indxstdppara[thisindxsampfull]]
            raise Exception('')
        gdatmodi.thisindxproptype = 0
    
    if gdat.booldiagmode and gdat.probspmr == 0 and gdatmodi.thisindxproptype > 2:
        print 'gdat.probtran'
        print gdat.probtran
        print 'gdat.probbrde'
        print gdat.probbrde
        raise Exception('')

    if gdat.verbtype > 1:
        print 'gdatmodi.thisindxproptype'
        print gdatmodi.thisindxproptype

    if numbtrap > 0:
        boolelemsbrtdfncanyy = getattr(gdat, strgmodl + 'boolelemsbrtdfncanyy')
        boolelemsbrtdfnc = getattr(gdat, strgmodl + 'boolelemsbrtdfnc')
        boolelemdeflsubh = getattr(gdat, strgmodl + 'boolelemdeflsubh')
        boolelemsbrtextsbgrd = getattr(gdat, strgmodl + 'boolelemsbrtextsbgrd')
        
    if gdatmodi.thisindxproptype == 0:
        nextsampunit = np.copy(thissampunit)
        if numbtrap > 0:
            nextindxelemfull = thisindxelemfull
    if gdatmodi.thisindxproptype > 0:
        nextsampunit = np.copy(thissampunit)
        nextsamp = np.copy(thissamp)
        if numbtrap > 0:
            nextindxelemfull = deepcopy(thisindxelemfull)
     
    if gdatmodi.thisindxproptype == 0:
        
        ## proposal scale
        if False:
            # amplitude-dependent proposal scale
            for l in indxpopl:
                thiscompampl = thissamp[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][indxelemfull]]
                compampl = nextsamp[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][indxelemfull]]
                minmcompampl = getattr(gdat, 'minm' + gdat.fittnamefeatampl[l])
                thiscompunit = thissamp[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][indxelemfull]]
                compunit = nextsamp[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][indxelemfull]]
                if strgcomp == gdat.fittnamefeatampl[l]:
                    # temp -- this only works if compampl is powr distributed
                    gdatmodi.thisstdp = stdpcomp / (thiscompampl / minmcompampl)**2.
                    gdatmodi.thisstdv = stdpcomp / (compampl / minmcompampl)**2.
                    gdatmodi.thislrpp += np.sum(0.5 * (nextcompunit - thiscompunit)**2 * (1. / gdatmodi.thisstdv**2 - 1. / gdatmodi.thisstdv**2))
                else:
                    gdatmodi.thisstdp = stdpcomp / (np.minimum(thiscompampl, compampl) / minmcompampl)**0.5
        
        ## propose a step
        diffsampunit = np.random.normal(size=thisindxsampfull.size) * thisstdp
        nextsampunit[thisindxsampfull] = thissampunit[thisindxsampfull] + diffsampunit
        
        if gdat.booldiagmode:
            if (nextsampunit[numbpopl:] == 1).any():
                raise Exception('')

            if (nextsampunit[numbpopl:] == 0).any():
                raise Exception('')

            if not np.isfinite(nextsampunit).all():
                print 'thisindxsampfull'
                summgene(thisindxsampfull)
                print 'thisstdp'
                print thisstdp
                summgene(thisstdp)
                raise Exception('')

        indxsamplowr = np.where(nextsampunit[numbpopl:] < 0.)[0]
        if indxsamplowr.size > 0:
            nextsampunit[numbpopl+indxsamplowr] = abs(nextsampunit[numbpopl+indxsamplowr]) % 1.
        
        if gdat.booldiagmode:
            if (nextsampunit[numbpopl:] == 1).any():
                raise Exception('')

            if (nextsampunit[numbpopl:] == 0).any():
                raise Exception('')

        indxsampuppr = np.where(nextsampunit[numbpopl:] > 1.)[0]
        if indxsampuppr.size > 0:
            nextsampunit[numbpopl+indxsampuppr] = (nextsampunit[numbpopl+indxsampuppr] - 1.) % 1.
        
        if gdat.booldiagmode:
            if (nextsampunit[numbpopl:] == 1).any():
                raise Exception('')

            if (nextsampunit[numbpopl:] == 0).any():
                raise Exception('')

            if not np.isfinite(nextsampunit).all():
                raise Exception('')

        nextsamp = icdf_samp(gdat, strgmodl, nextsampunit, thisindxsampcomp)

    if gdatmodi.thisindxproptype > 0:
        gdatmodi.indxsamptran = []
        if gdatmodi.thisindxproptype == 1:
            gdatmodi.thisauxipara = rand(numbcomp[gdatmodi.indxpopltran])
        elif gdatmodi.thisindxproptype != 2:
            gdatmodi.thisauxipara = np.empty(numbcomp[gdatmodi.indxpopltran])
    
    if gdatmodi.thisindxproptype == 1 or gdatmodi.thisindxproptype == 3:
       
        # find an empty slot in the element list
        for u in range(maxmnumbelem[gdatmodi.indxpopltran]):
            if not u in gdatmodi.thisindxelemfull[gdatmodi.indxpopltran]:
                break
        gdatmodi.indxelemmodi = [u]
        gdatmodi.indxelemfullmodi = [thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]].astype(int)]
       
        # sample indices to add the new element
        gdatmodi.indxsampelemaddd = retr_indxsamppnts(indxsampcompinit, numbtrappoplcuml, numbcomp, indxcomp, \
                                                                                                       gdatmodi.indxpopltran, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.append(gdatmodi.indxsampelemaddd)
        nextindxelemfull[gdatmodi.indxpopltran].append(gdatmodi.indxelemmodi[0])
    if gdatmodi.thisindxproptype == 1:
        
        # sample auxiliary variables
        nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.thisauxipara
    
    # death
    if gdatmodi.thisindxproptype == 2:
        
        # occupied element index to be killed
        if thisindxelem is None:
            dethindxindxelem = choice(np.arange(thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]], dtype=int))
        else:
            dethindxindxelem = thisindxelem

        # element index to be killed
        gdatmodi.indxelemmodi = []
        gdatmodi.indxelemfullmodi = []
        if gdat.verbtype > 1:
            print 'dethindxindxelem'
            print dethindxindxelem

        gdatmodi.indxelemmodi.append(thisindxelemfull[gdatmodi.indxpopltran][dethindxindxelem])
        gdatmodi.indxelemfullmodi.append(dethindxindxelem)
        # parameter indices to be killed
        indxsampelemdeth = retr_indxsamppnts(indxsampcompinit, numbtrappoplcuml, numbcomp, indxcomp, \
                                                                                    gdatmodi.indxpopltran, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.append(indxsampelemdeth)
        
        gdatmodi.thisauxipara = thissamp[indxsampelemdeth]

    if gdatmodi.thisindxproptype > 2:
        gdatmodi.comppare = np.empty(numbcomp[gdatmodi.indxpopltran])
        gdatmodi.compfrst = np.empty(numbcomp[gdatmodi.indxpopltran])
        gdatmodi.compseco = np.empty(numbcomp[gdatmodi.indxpopltran])
    
    # split
    if gdatmodi.thisindxproptype == 3:
        
        # find the probability of splitting elements
        probsplt = retr_probspmr(gdat, gdatmodi, thissamp, thisindxsampcomp, gdatmodi.indxpopltran, boolmerg=False)
        gdatmodi.indxelemfullsplt = choice(np.arange(thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]], dtype=int), p=probsplt)
        gdatmodi.indxelemsplt = thisindxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullsplt]
        gdatmodi.indxelemfullmodi.insert(0, gdatmodi.indxelemfullsplt)
        gdatmodi.indxelemmodi.insert(0, gdatmodi.indxelemsplt)

        # sample indices for the first element
        gdatmodi.indxsampelemfrst = retr_indxsamppnts(indxsampcompinit, numbtrappoplcuml, numbcomp, indxcomp, \
                                                                    l, gdatmodi.indxelemmodi[0])
        gdatmodi.indxsamptran.insert(0, gdatmodi.indxsampelemfrst)
        
        # sample indices for the second element
        gdatmodi.indxsampseco = gdatmodi.indxsampelemaddd
        
        # take the parent element parameters
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            gdatmodi.comppare[k] = np.copy(thissamp[thisindxsampcomp[strgcomp][gdatmodi.indxpopltran][gdatmodi.indxelemfullmodi[0]]])
        
        # draw the auxiliary parameters
        for g, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            if boolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.thisauxipara[g] = randn() * gdat.radispmr
            elif g == indxcompampl[gdatmodi.indxpopltran]:
                gdatmodi.thisauxipara[g] = rand()
            else:
                gdatmodi.thisauxipara[g] = icdf_trap(gdat, strgmodl, np.random.rand(), thissamp, listscalcomp[gdatmodi.indxpopltran][g], \
                                                                            liststrgcomp[gdatmodi.indxpopltran][g], l)

        # determine the new parameters
        if elemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.thisauxipara[1]) * gdatmodi.thisauxipara[0]
        else:
            gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.thisauxipara[2]) * gdatmodi.thisauxipara[0]
            gdatmodi.compfrst[1] = gdatmodi.comppare[1] + (1. - gdatmodi.thisauxipara[2]) * gdatmodi.thisauxipara[1]
        gdatmodi.compfrst[indxcompampl[gdatmodi.indxpopltran]] = gdatmodi.thisauxipara[indxcompampl[gdatmodi.indxpopltran]] * \
                                                                                                        gdatmodi.comppare[indxcompampl[gdatmodi.indxpopltran]]
        if elemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.thisauxipara[1] * gdatmodi.thisauxipara[0]
        else:
            gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.thisauxipara[2] * gdatmodi.thisauxipara[0]
            gdatmodi.compseco[1] = gdatmodi.comppare[1] - gdatmodi.thisauxipara[2] * gdatmodi.thisauxipara[1]
        gdatmodi.compseco[indxcompampl[gdatmodi.indxpopltran]] = (1. - gdatmodi.thisauxipara[indxcompampl[gdatmodi.indxpopltran]]) * \
                                                                                                        gdatmodi.comppare[indxcompampl[gdatmodi.indxpopltran]]
        for g in range(numbcomp[gdatmodi.indxpopltran]):
            if not boolcompposi[gdatmodi.indxpopltran][g] and g != indxcompampl[gdatmodi.indxpopltran]:
                gdatmodi.compfrst[g] = gdatmodi.comppare[g]
                gdatmodi.compseco[g] = gdatmodi.thisauxipara[g]
       
        # place the new parameters into the sample vector
        nextsamp[gdatmodi.indxsamptran[0]] = cdfn_trap(gdat, gdatmodi, gdatmodi.compfrst, gdatmodi.indxpopltran)
        nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.compfrst
        nextsamp[gdatmodi.indxsamptran[1]] = cdfn_trap(gdat, gdatmodi, gdatmodi.compseco, gdatmodi.indxpopltran)
        nextsamp[gdatmodi.indxsamptran[1]] = gdatmodi.compseco
        
        # check for prior boundaries
        if elemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            if np.fabs(gdatmodi.compfrst[0]) > gdat.maxmelin or np.fabs(gdatmodi.compseco[0]) > gdat.maxmelin:
                gdatmodi.thisboolpropfilt = False
        else:
            maxmlgal = getattr(gdat, strgmodl + 'maxmlgal')
            maxmbgal = getattr(gdat, strgmodl + 'maxmbgal')
            if np.fabs(gdatmodi.compfrst[0]) > maxmlgal or np.fabs(gdatmodi.compseco[0]) > maxmlgal or \
                                                                    np.fabs(gdatmodi.compfrst[1]) > maxmbgal or np.fabs(gdatmodi.compseco[1]) > maxmbgal:
                gdatmodi.thisboolpropfilt = False
        if gdatmodi.compfrst[indxcompampl[gdatmodi.indxpopltran]] < getattr(gdat, strgmodl + 'minm' + gdat.fittnamefeatampl[gdatmodi.indxpopltran]) or \
                  gdatmodi.compseco[indxcompampl[gdatmodi.indxpopltran]] < getattr(gdat, strgmodl + 'minm' + gdat.fittnamefeatampl[gdatmodi.indxpopltran]):
            gdatmodi.thisboolpropfilt = False
        if gdat.verbtype > 1:
            if not gdatmodi.thisboolpropfilt:
                print 'Rejecting the proposal due to a split that falls out of the prior...'
    
    if gdatmodi.thisindxproptype == 4:
        
        # determine the index of the primary element to be merged (in the full element list)
        gdatmodi.indxelemfullmergfrst = choice(np.arange(len(thisindxelemfull[gdatmodi.indxpopltran])))

        ## first element index to be merged
        gdatmodi.mergindxelemfrst = thisindxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullmergfrst]
         
        # find the probability of merging this element with the others 
        probmerg = retr_probspmr(gdat, gdatmodi, thissamp, thisindxsampcomp, gdatmodi.indxpopltran, \
                            indxelemfullexcl=gdatmodi.indxelemfullmergfrst, elemtype=elemtype)
        
        indxelemfulltemp = np.arange(len(thisindxelemfull[gdatmodi.indxpopltran]))
        if gdat.booldiagmode:
            if indxelemfulltemp.size < 2:
                print 'indxelemfulltemp'
                print indxelemfulltemp
                raise Exception('')
        gdatmodi.indxelemfullmergseco = choice(np.setdiff1d(indxelemfulltemp, np.array([gdatmodi.indxelemfullmergfrst])), p=probmerg)
        gdatmodi.indxelemfullmodi = np.sort(np.array([gdatmodi.indxelemfullmergfrst, gdatmodi.indxelemfullmergseco]))
        
        # parameters of the first element to be merged
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            ## first
            gdatmodi.compfrst[k] = thissamp[thisindxsampcomp[strgcomp][gdatmodi.indxpopltran][gdatmodi.indxelemfullmodi[0]]]
        
        # determine indices of the modified elements in the sample vector
        ## first element
        # temp -- this would not work for multiple populations !
        gdatmodi.indxsampelemfrst = retr_indxsamppnts(indxsampcompinit, numbtrappoplcuml, numbcomp, indxcomp, \
                                                                    l, gdatmodi.mergindxelemfrst)
        gdatmodi.indxsamptran.append(gdatmodi.indxsampelemfrst)

        ## second element index to be merged
        gdatmodi.mergindxelemseco = thisindxelemfull[gdatmodi.indxpopltran][gdatmodi.indxelemfullmergseco]
       
        ## second element
        gdatmodi.indxsampelemseco = retr_indxsamppnts(indxsampcompinit, numbtrappoplcuml, numbcomp, indxcomp, \
                                                                    l, gdatmodi.mergindxelemseco)
        gdatmodi.indxsamptran.append(gdatmodi.indxsampelemseco)
        
        # parameters of the elements to be merged
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            ## second
            gdatmodi.compseco[k] = thissamp[thisindxsampcomp[strgcomp][gdatmodi.indxpopltran][gdatmodi.indxelemfullmodi[1]]]

        # indices of the element to be merged
        gdatmodi.indxelemmodi = [gdatmodi.mergindxelemfrst, gdatmodi.mergindxelemseco]

        # auxiliary parameters
        if elemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.thisauxipara[0] = gdatmodi.compseco[0] - gdatmodi.compfrst[0]
        else:
            gdatmodi.thisauxipara[0] = gdatmodi.compseco[0] - gdatmodi.compfrst[0]
            gdatmodi.thisauxipara[1] = gdatmodi.compseco[1] - gdatmodi.compfrst[1]
        gdatmodi.thisauxipara[indxcompampl[gdatmodi.indxpopltran]] = gdatmodi.compfrst[indxcompampl[gdatmodi.indxpopltran]] / \
                                        (gdatmodi.compfrst[indxcompampl[gdatmodi.indxpopltran]] + gdatmodi.compseco[indxcompampl[gdatmodi.indxpopltran]]) 
        for g, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            if not boolcompposi[gdatmodi.indxpopltran][g] and g != indxcompampl[gdatmodi.indxpopltran]:
                gdatmodi.thisauxipara[g] = gdatmodi.compseco[g]

        # merged element
        gdatmodi.comppare[indxcompampl[gdatmodi.indxpopltran]] = gdatmodi.compfrst[indxcompampl[gdatmodi.indxpopltran]] + \
                                                                                                gdatmodi.compseco[indxcompampl[gdatmodi.indxpopltran]]
        if gdatmodi.comppare[indxcompampl[gdatmodi.indxpopltran]] > getattr(gdat, 'maxm' + gdat.fittnamefeatampl[gdatmodi.indxpopltran]):
            gdatmodi.thisboolpropfilt = False
            if gdat.verbtype > 1:
                print 'Proposal rejected due to falling outside the prior.'
            return

        if elemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.comppare[0] = gdatmodi.compfrst[0] + (1. - gdatmodi.thisauxipara[1]) * (gdatmodi.compseco[0] - gdatmodi.compfrst[0])
        else:
            gdatmodi.comppare[0] = gdatmodi.compfrst[0] + (1. - gdatmodi.thisauxipara[2]) * (gdatmodi.compseco[0] - gdatmodi.compfrst[0])
            gdatmodi.comppare[1] = gdatmodi.compfrst[1] + (1. - gdatmodi.thisauxipara[2]) * (gdatmodi.compseco[1] - gdatmodi.compfrst[1])
        for g, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            if boolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.comppare[g] = gdatmodi.compfrst[g] + (1. - gdatmodi.thisauxipara[indxcompampl[gdatmodi.indxpopltran]]) * \
                                                                                            (gdatmodi.compseco[g] - gdatmodi.compfrst[g])
            elif g == indxcompampl[gdatmodi.indxpopltran]:
                gdatmodi.comppare[g] = gdatmodi.compfrst[g] + gdatmodi.compseco[g]
            else:
                gdatmodi.comppare[g] = gdatmodi.compfrst[g]

        nextsamp[gdatmodi.indxsamptran[0]] = cdfn_trap(gdat, gdatmodi, gdatmodi.comppare, gdatmodi.indxpopltran)
        nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.comppare

        # calculate the proposed list of pairs
        if gdat.verbtype > 1:
            print 'mergindxfrst: ', gdatmodi.mergindxelemfrst
            print 'gdatmodi.indxelemfullmergfrst: ', gdatmodi.indxelemfullmergfrst
            print 'mergindxseco: ', gdatmodi.mergindxelemseco
            print 'gdatmodi.indxelemfullmergseco: ', gdatmodi.indxelemfullmergseco
            print 'indxsampelemfrst: ', gdatmodi.indxsampelemfrst
            print 'indxsampelemseco: ', gdatmodi.indxsampelemseco

    if gdat.verbtype > 1 and (gdatmodi.thisindxproptype == 3 or gdatmodi.thisboolpropfilt and gdatmodi.thisindxproptype == 4):
        
        if elemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            print 'elinfrst: ', gdatmodi.compfrst[0]
            print 'amplfrst: ', gdatmodi.compfrst[1]
            print 'elinseco: ', gdatmodi.compseco[0]
            print 'amplseco: ', gdatmodi.compseco[1]
            print 'elinpare: ', gdatmodi.comppare[0]
            print 'fluxpare: ', gdatmodi.comppare[1]
            print 'auxipara[0][0]: ', gdatmodi.thisauxipara[0]
            print 'auxipara[0][1]: ', gdatmodi.thisauxipara[1]
        else:
            print 'lgalfrst: ', gdat.anglfact * gdatmodi.compfrst[0]
            print 'bgalfrst: ', gdat.anglfact * gdatmodi.compfrst[1]
            print 'amplfrst: ', gdatmodi.compfrst[2]
            print 'lgalseco: ', gdat.anglfact * gdatmodi.compseco[0]
            print 'bgalseco: ', gdat.anglfact * gdatmodi.compseco[1]
            print 'amplseco: ', gdatmodi.compseco[2]
            print 'lgalpare: ', gdat.anglfact * gdatmodi.comppare[0]
            print 'bgalpare: ', gdat.anglfact * gdatmodi.comppare[1]
            print 'fluxpare: ', gdatmodi.comppare[2]
            print 'auxipara[0][0]: ', gdat.anglfact * gdatmodi.thisauxipara[0]
            print 'auxipara[0][1]: ', gdat.anglfact * gdatmodi.thisauxipara[1]
            print 'auxipara[0][2]: ', gdatmodi.thisauxipara[2]
                
    if numbtrap > 0 and gdatmodi.thisindxproptype > 0 and gdatmodi.thisboolpropfilt:
        # change the number of elements
        if gdatmodi.thisindxproptype == 1 or gdatmodi.thisindxproptype == 3:
            nextsamp[indxfixpnumbelem[gdatmodi.indxpopltran]] = thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]] + 1
        if gdatmodi.thisindxproptype == 2 or gdatmodi.thisindxproptype == 4:
            nextsamp[indxfixpnumbelem[gdatmodi.indxpopltran]] = thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]] - 1
        nextsampunit[indxfixpnumbelem[gdatmodi.indxpopltran]] = nextsamp[indxfixpnumbelem[gdatmodi.indxpopltran]]
        
        # remove the element from the occupied element list
        if (gdatmodi.thisindxproptype == 2 or gdatmodi.thisindxproptype == 4):
            for a, indxelem in enumerate(gdatmodi.indxelemmodi):
                if a == 0 and gdatmodi.thisindxproptype == 2 or a == 1 and gdatmodi.thisindxproptype == 4:
                    nextindxelemfull[gdatmodi.indxpopltran].remove(indxelem)
    
    if gdatmodi.thisindxproptype == 0:
        gdatmodi.indxsampmodi = thisindxsampfull
    else:
        if gdatmodi.thisindxproptype == 1:
            gdatmodi.indxsampmodi = np.concatenate((np.array([indxfixpnumbelem[gdatmodi.indxpopltran]]), gdatmodi.indxsamptran[0]))
        if gdatmodi.thisindxproptype == 2:
            gdatmodi.indxsampmodi = [indxfixpnumbelem[gdatmodi.indxpopltran]]
        if gdatmodi.thisindxproptype == 3:
            gdatmodi.indxsampmodi = np.concatenate((np.array([indxfixpnumbelem[gdatmodi.indxpopltran]]), \
                                                                            gdatmodi.indxsamptran[0], gdatmodi.indxsamptran[1]))
        if gdatmodi.thisindxproptype == 4:
            gdatmodi.indxsampmodi = np.concatenate((np.array([indxfixpnumbelem[gdatmodi.indxpopltran]]), gdatmodi.indxsamptran[0]))
    
    if numbtrap > 0:
        if gdatmodi.thisindxproptype == 0:
            indxsampcomp = thisindxsampcomp
        else:
            indxsampcomp = retr_indxsampcomp(gdat, nextindxelemfull, strgmodl)
    if gdat.verbtype > 1:
        print 'gdatmodi.indxsampmodi'
        print gdatmodi.indxsampmodi
        if numbtrap > 0:
            print 'thisindxelemfull'
            print thisindxelemfull
            print 'thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]].astype(int)'
            print thissamp[indxfixpnumbelem[gdatmodi.indxpopltran]].astype(int)
            if gdatmodi.thisindxproptype > 0:
                print 'gdatmodi.indxelemmodi'
                print gdatmodi.indxelemmodi
                print 'gdatmodi.indxelemfullmodi'
                print gdatmodi.indxelemfullmodi
                print 'gdatmodi.thisboolpropfilt'
                print gdatmodi.thisboolpropfilt
            print 'indxsampcomp'
            print indxsampcomp
    
    if gdatmodi.thisindxproptype == 1:
        for g, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpopltran]):
            nextsamp[gdatmodi.indxsamptran[0][g]] = icdf_trap(gdat, strgmodl, \
                gdatmodi.thisauxipara[g], thissamp, listscalcomp[gdatmodi.indxpopltran][g], liststrgcomp[gdatmodi.indxpopltran][g], gdatmodi.indxpopltran)

    setattr(gdatobjt, strgpfixnext + 'sampunit', nextsampunit)
    setattr(gdatobjt, strgpfixnext + 'samp', nextsamp)
    if numbtrap > 0:
        setattr(gdatobjt, strgpfixnext + 'indxsampcomp', indxsampcomp)
        setattr(gdatobjt, strgpfixnext + 'indxelemfull', nextindxelemfull)
        
    if gdat.booldiagmode:
        if numbtrap > 0:
            for l in indxpopl:
                if thissampunit[indxfixpnumbelem[l]] != round(thissampunit[indxfixpnumbelem[l]]):
                    raise Exception('')
                if thissamp[indxfixpnumbelem[l]] != round(thissamp[indxfixpnumbelem[l]]):
                    raise Exception('')
                if nextsampunit[indxfixpnumbelem[l]] != round(nextsampunit[indxfixpnumbelem[l]]):
                    raise Exception('')
                if nextsamp[indxfixpnumbelem[l]] != round(nextsamp[indxfixpnumbelem[l]]):
                    raise Exception('')

        if strgmodl == 'fitt':
            diffsamp = abs(nextsamp - thissamp)
            #size = np.where(((thissamp == 0.) & (diffsamp > 0.)) | ((thissamp != 0.) & (diffsamp / thissamp > 0)))[0].size
            size = np.where(diffsamp != 0.)[0].size
            if gdatmodi.thisindxproptype == 1:
                if size - 1 != numbcomp[gdatmodi.indxpopltran]:
                    print 'Changing elsenp.where!!'
                    print 'gdatmodi.thisindxproptype'
                    print gdatmodi.thisindxproptype
                    show_samp(gdat, gdatmodi, 'this', 'fitt')
                    raise Exception('')
    

def calc_probprop(gdat, gdatmodi):

    # calculate the factor to multiply the acceptance rate, i.e., 
    ## probability of the auxiliary parameters,
    if gdatmodi.thisindxproptype == 0:
        gdatmodi.thislpau = 0.
    elif gdatmodi.thisindxproptype == 1 or gdatmodi.thisindxproptype == 2:
        gdatmodi.thislpau = gdatmodi.nextlpritotl - gdatmodi.thislpritotl
        lpautemp = 0.5 * gdat.priofactdoff * gdat.fittnumbcomp[gdatmodi.indxpopltran]
        if gdatmodi.thisindxproptype == 1:
            gdatmodi.thislpau += lpautemp
        if gdatmodi.thisindxproptype == 2:
            gdatmodi.thislpau -= lpautemp
    elif gdatmodi.thisindxproptype == 3 or gdatmodi.thisindxproptype == 4:
        gdatmodi.thislpau = 0.
        liststrgpdfnprio = getattr(gdat, 'fittliststrgpdfnprio')
        dictelemtemp = [dict()]
        for g, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpopltran]):
            if gdat.fittboolcompposi[gdatmodi.indxpopltran][g]:
                gdatmodi.thislpau += -0.5 * np.log(2. * np.pi * gdat.radispmr**2) - 0.5 * (gdatmodi.thisauxipara[g] / gdat.radispmr)**2
            elif g != gdat.fittindxcompampl[gdatmodi.indxpopltran]:
                dictelemtemp[0][strgcomp] = gdatmodi.thisauxipara[g]
                gdatmodi.thislpau += retr_lprielem(gdat, 'fitt', gdatmodi.indxpopltran, g, \
                                            gdat.fittliststrgcomp[gdatmodi.indxpopltran][g], gdat.fittliststrgpdfnprio[gdatmodi.indxpopltran][g], \
                                            gdatmodi.thissamp, dictelemtemp, [1])
        if gdatmodi.thisindxproptype == 4:
            gdatmodi.thislpau *= -1.

    if gdatmodi.thisindxproptype > 2 and gdatmodi.thisboolpropfilt:
        ## the ratio of the probability of the reverse and forward proposals, and
        if gdatmodi.thisindxproptype == 3:
            gdatmodi.thisprobrevefrst = retr_probspmr(gdat, gdatmodi, gdatmodi.nextsamp, gdatmodi.nextindxsampcomp, gdatmodi.indxpopltran, \
                                                                               indxelemfullexcl=-1, elemtype=gdat.fittelemtype, \
                                                                               indxelemfulleval=gdatmodi.indxelemfullmodi[0])
            gdatmodi.thisprobreveseco = retr_probspmr(gdat, gdatmodi, gdatmodi.nextsamp, gdatmodi.nextindxsampcomp, gdatmodi.indxpopltran, \
                                             indxelemfullexcl=gdatmodi.indxelemfullmodi[0], elemtype=gdat.fittelemtype, indxelemfulleval=-1)
            gdatmodi.thislrpp = np.log(gdatmodi.thisprobrevefrst + gdatmodi.thisprobreveseco) \
                                    - np.log(retr_probspmr(gdat, gdatmodi, gdatmodi.nextsamp, gdatmodi.nextindxsampcomp, gdatmodi.indxpopltran, boolmerg=False, \
                                      indxelemfulleval=gdatmodi.indxelemfullmodi[0]))

        else:
            gdatmodi.thisprobfrwdfrst = retr_probspmr(gdat, gdatmodi, gdatmodi.thissamp, gdatmodi.thisindxsampcomp, gdatmodi.indxpopltran, \
                                                                            indxelemfullexcl=gdatmodi.indxelemfullmodi[1], \
                                                                            elemtype=gdat.fittelemtype, indxelemfulleval=gdatmodi.indxelemfullmodi[0])
            gdatmodi.thisprobfrwdseco = retr_probspmr(gdat, gdatmodi, gdatmodi.thissamp, gdatmodi.thisindxsampcomp, gdatmodi.indxpopltran, \
                                                                            indxelemfullexcl=gdatmodi.indxelemfullmodi[0], \
                                                                            elemtype=gdat.fittelemtype, indxelemfulleval=gdatmodi.indxelemfullmodi[1])
            
            gdatmodi.thislrpp = -np.log(gdatmodi.thisprobfrwdfrst + gdatmodi.thisprobfrwdseco) \
                                    - np.log(retr_probspmr(gdat, gdatmodi, gdatmodi.nextsamp, gdatmodi.nextindxsampcomp, gdatmodi.indxpopltran, boolmerg=False, \
                                          indxelemfulleval=gdatmodi.indxelemfullmodi[0]))
       
        ## Jacobian
        if gdat.fittelemtype[gdatmodi.indxpopltran].startswith('lghtline'):
            gdatmodi.thisljcb = np.log(gdatmodi.comppare[1])
        else:
            gdatmodi.thisljcb = np.log(gdatmodi.comppare[2])
        if gdatmodi.thisindxproptype == 4:
            gdatmodi.thisljcb *= -1.
        
    else:
        gdatmodi.thisljcb = 0.
        gdatmodi.thislrpp = 0.
    
    for l in gdat.fittindxpopl:
        if gdatmodi.thisindxproptype > 0:
            setattr(gdatmodi, 'auxiparapop%d' % l, gdatmodi.thisauxipara)


def retr_indxsampcomp(gdat, indxelemfull, strgmodl):

    ## element parameters
    #strgpfix = retr_strgpfix(strgstat, strgmodl)
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        #indxelemfull = list(getattr(gdatobjt, strgpfix + 'indxelemfull'))
        numbtrappoplcuml = getattr(gdat, strgmodl + 'numbtrappoplcuml')
        elemtype = getattr(gdat, strgmodl + 'elemtype')
        numbcomp = getattr(gdat, strgmodl + 'numbcomp')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        boolelemlght = getattr(gdat, strgmodl + 'boolelemlght')
        boolelemlghtspat = getattr(gdat, strgmodl + 'boolelemlghtspat')
        boolelemsbrtextsbgrd = getattr(gdat, strgmodl + 'boolelemsbrtextsbgrd')
        spatdisttype = getattr(gdat, strgmodl + 'spatdisttype')
        indxcomp = getattr(gdat, strgmodl + 'indxcomp')
        spectype = getattr(gdat, strgmodl + 'spectype')
        indxsampcompinit = getattr(gdat, strgmodl + 'indxsampcompinit')
        liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
        
        indxsampcomp = dict()
        for strgcomp in liststrgcomptotl:
            indxsampcomp[strgcomp] = [[] for l in indxpopl]

        indxsampcomp['comp'] = [[] for l in indxpopl]
        for l in indxpopl:
            indxsamptemp = indxsampcompinit + numbtrappoplcuml[l] + np.array(indxelemfull[l], dtype=int) * numbcomp[l]
            cntr = tdpy.util.cntr()
            # position
            if elemtype[l][:8] == 'lghtline':
                indxsampcomp['elin'][l] = indxsamptemp + cntr.incr()
            elif spatdisttype[l] == 'gangexpo':
                indxsampcomp['gang'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['aang'][l] = indxsamptemp + cntr.incr()
            elif spatdisttype[l] == 'glc3':
                indxsampcomp['dglc'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['thet'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['phii'][l] = indxsamptemp + cntr.incr()
            elif spatdisttype[l] == 'los3':
                indxsampcomp['lgal'][l] = indxsamptemp + cntr.incr() 
                indxsampcomp['bgal'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['dlos'][l] = indxsamptemp + cntr.incr()
            else:
                indxsampcomp['lgal'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['bgal'][l] = indxsamptemp + cntr.incr()
            
            # amplitude
            if elemtype[l] == 'lghtpntspuls':
                indxsampcomp['per0'][l] = indxsamptemp + cntr.incr()
            elif elemtype[l] == 'lghtpntsagnntrue':
                indxsampcomp['lum0'][l] = indxsamptemp + cntr.incr()
            elif boolelemlght[l]:
                indxsampcomp['flux'][l] = indxsamptemp + cntr.incr()
            elif elemtype[l] == 'lens':
                indxsampcomp['defs'][l] = indxsamptemp + cntr.incr()
            elif elemtype[l].startswith('clus'):
                indxsampcomp['nobj'][l] = indxsamptemp + cntr.incr()
            
            # shape
            if boolelemsbrtextsbgrd[l] or elemtype[l] == 'clusvari':
                indxsampcomp['gwdt'][l] = indxsamptemp + cntr.incr()
            if elemtype[l] == 'lens':
                if gdat.variasca:
                    indxsampcomp['asca'][l] = indxsamptemp + cntr.incr()
                if gdat.variacut:
                    indxsampcomp['acut'][l] = indxsamptemp + cntr.incr()
            if elemtype[l] == 'lghtlinevoig':
                indxsampcomp['sigm'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['gamm'][l] = indxsamptemp + cntr.incr()

            # others
            if elemtype[l] == 'lghtpntspuls':
                indxsampcomp['magf'][l] = indxsamptemp + cntr.incr()
                indxsampcomp['geff'][l] = indxsamptemp + cntr.incr()
            elif elemtype[l] == 'lghtpntsagnntrue':
                indxsampcomp['dlos'][l] = indxsamptemp + cntr.incr()
            
            # spectra
            if boolelemlghtspat[l]:
                if gdat.numbener > 1:
                    if spectype[l] == 'colr':
                        for i in gdat.indxener:
                            if i == 0:
                                continue
                            indxsampcomp['sindcolr%04d' % i][l] = indxsamptemp + cntr.incr()
                    else:
                        indxsampcomp['sind'][l] = indxsamptemp + cntr.incr()
                        if spectype[l] == 'curv':
                            indxsampcomp['curv'][l] = indxsamptemp + cntr.incr()
                        if spectype[l] == 'expc':
                            indxsampcomp['expc'][l] = indxsamptemp + cntr.incr()
            
            indxsampcomp['comp'][l] = np.repeat(indxsamptemp, numbcomp[l]) + np.tile(indxcomp[l], len(indxelemfull[l]))
        
        if gdat.booldiagmode:
            numbpara = getattr(gdat, strgmodl + 'numbpara')
            for l in indxpopl:
                if len(indxsampcomp['comp'][l]) > 0:
                    if np.amax(indxsampcomp['comp'][l]) > numbpara:
                        print 'strgmodl'
                        print strgmodl
                        print 'l'
                        print l
                        print 'indxsampcomp[comp][l]'
                        summgene(indxsampcomp['comp'][l])
                        print 'numbpara'
                        print numbpara
                        #print 'samp'
                        #for k in range(samp.size):
                        #    print samp[k]
                        raise Exception('Sample indices of the elements are bad.') 
        
    else:
        indxsampcomp = None
    
    return indxsampcomp
    

def retr_probspmr(gdat, gdatmodi, samp, indxsampcomp, indxpopltran, indxelemfullexcl=None, elemtype=None, indxelemfulleval=None, boolmerg=True):
    
    if boolmerg:
        if elemtype[indxpopltran].startswith('lghtline'):
            elin = samp[indxsampcomp['elin'][indxpopltran]][indxelemfullexcl]
            elinstat = samp[indxsampcomp['elin'][indxpopltran]]
            probmerg = np.exp(-0.5 * ((elin - elinstat) / gdat.radispmr)**2)
        else:
            lgal = samp[indxsampcomp['lgal'][indxpopltran]][indxelemfullexcl]
            bgal = samp[indxsampcomp['bgal'][indxpopltran]][indxelemfullexcl]
            lgalstat = samp[indxsampcomp['lgal'][indxpopltran]]
            bgalstat = samp[indxsampcomp['bgal'][indxpopltran]]
            probmerg = np.exp(-0.5 * (((lgal - lgalstat) / gdat.radispmr)**2 + ((bgal - bgalstat) / gdat.radispmr)**2))
        
        if (probmerg == 0.).all():
            gdatmodi.thisboolpropfilt = False
            print 'Merhe probability is zero.'
            print 'Turning off gdatmodi.thisboolpropfilt'
            probmerg[0] = 1.
            print 'probmerg'
            print probmerg
            print

    else:
        amplstat = samp[indxsampcomp[gdat.fittnamefeatampl[indxpopltran]][indxpopltran]]
        probmerg = amplstat / np.sum(amplstat)
    
    if indxelemfulleval is None:
        if indxelemfullexcl is not None:
            probmerg = np.concatenate((probmerg[:indxelemfullexcl], probmerg[indxelemfullexcl+1:]))
            probmerg /= np.sum(probmerg)
    else:
        if indxelemfullexcl is not None:
            if indxelemfullexcl == -1:
                probmerg = probmerg[:-1]
            else:
                probmerg = np.concatenate((probmerg[:indxelemfullexcl], probmerg[indxelemfullexcl+1:]))
            probmerg /= np.sum(probmerg)
        if indxelemfullexcl < indxelemfulleval and indxelemfulleval != -1:
            indxelemfulleval -= 1
        probmerg = probmerg[indxelemfulleval]

    if not np.isfinite(probmerg).all():
        gdatmodi.thisboolpropfilt = False
        probmerg = np.zeros_like(probmerg)
        probmerg[0] = 1.

    if gdat.booldiagmode:
        if not np.isfinite(probmerg).all():
            print 'Merge probability is infinite.'
            print 'probmerg'
            print probmerg
            print 'boolmerg'
            print boolmerg
            print 'indxelemfullexcl'
            print indxelemfullexcl
            print 'indxelemfulleval'
            print indxelemfulleval
            print
            #raise Exception('Merge probability is infinite.')

    return probmerg

    
def retr_indxsamppnts(indxsampcompinit, numbtrappoplcuml, numbcomp, indxcomp, l, u):

    indxsamppnts = indxsampcompinit + numbtrappoplcuml[l] + u * numbcomp[l] + indxcomp[l]

    return indxsamppnts


def gang_detr():

    gang, aang, lgal, bgal = sympy.symbols('gang aang lgal bgal')

    AB = sympy.matrices.Matrix([[a1*b1,a1*b2,a1*b3],[a2*b1,a2*b2,a2*b3],[a3*b1,a3*b2,a3*b3]])
    print AB.det()


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, strgmodl):

    numbpsfpform = getattr(gdat, strgmodl + 'numbpsfpform')
    numbpsfptotl = getattr(gdat, strgmodl + 'numbpsfptotl')
    
    indxpsfpinit = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * np.arcsin(np.sqrt(2. - 2. * np.cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * np.arcsin(np.sqrt(2. - 2. * np.cos(gdat.binsangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        scalangl = thisangl[None, :, None]
    
    if psfntype == 'singgaus':
        sigc = psfp[indxpsfpinit]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
    
    elif psfntype == 'singking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(scalangl, sigc, gamc)
    
    elif psfntype == 'doubking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        sigt = psfp[indxpsfpinit+2]
        gamt = psfp[indxpsfpinit+3]
        frac = psfp[indxpsfpinit+4]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        frac = frac[:, None, :]
        
        psfn = retr_doubking(scalangl, frac, sigc, gamc, sigt, gamt)
        if gdat.exprtype == 'ferm':
            psfnnorm = retr_doubking(scalanglnorm, frac, sigc, gamc, sigt, gamt)
    
    # normalize the PSF
    if gdat.exprtype == 'ferm':
        fact = 2. * np.pi * np.trapz(psfnnorm * np.sin(gdat.binsangl[None, :, None]), gdat.binsangl, axis=1)[:, None, :]
        psfn /= fact

    return psfn


def retr_axis(gdat, strgvarb, minm=None, maxm=None, numb=None, bins=None, scal='self', invr=False, strginit=''):
    
    if bins is None:
        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            binsscal = np.linspace(minm, maxm, numb + 1)
        if scal == 'logt':
            binsscal = np.linspace(np.log10(minm), np.log10(maxm), numb + 1)
        if scal == 'asnh':
            binsscal = np.linspace(np.arcsinh(minm), np.arcsinh(maxm), numb + 1)
        
        if invr:
            binsscal = binsscal[::-1]
        
        meanvarbscal = (binsscal[1:] + binsscal[:-1]) / 2.
        
        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            meanvarb = meanvarbscal
            bins = binsscal
        if scal == 'logt':
            meanvarb = 10**meanvarbscal
            bins = 10**binsscal
        if scal == 'asnh':
            meanvarb = np.sinh(meanvarbscal)
            bins = np.sinh(binsscal)
    else:
        numb = bins.size - 1

    indx = np.arange(numb)
    delt = np.diff(bins) 
    limt = np.array([np.amin(bins), np.amax(bins)]) 

    setattr(gdat, strginit + 'limt' + strgvarb, limt)
    setattr(gdat, strginit + 'bins' + strgvarb, bins)
    setattr(gdat, strginit + 'mean' + strgvarb, meanvarb)
    setattr(gdat, strginit + 'delt' + strgvarb, delt)
    setattr(gdat, strginit + 'numb' + strgvarb, numb)
    setattr(gdat, strginit + 'indx' + strgvarb, indx)


def retr_unit(lgal, bgal):

    xdat = np.cos(bgal) * np.cos(lgal)
    ydat = -np.cos(bgal) * np.sin(lgal)
    zaxi = np.sin(bgal)

    return xdat, ydat, zaxi


def retr_psec(gdat, conv):

    # temp
    conv = conv.reshape((gdat.numbsidecart, gdat.numbsidecart))
    psec = (abs(scipy.fftpack.fft2(conv))**2)[:gdat.numbsidecart/2, :gdat.numbsidecart/2] * 1e-3
    psec = psec.flatten()

    return psec
   

def retr_psecodim(gdat, psec):
    
    psec = psec.reshape((gdat.numbsidecart / 2, gdat.numbsidecart / 2))
    psecodim = np.zeros(gdat.numbsidecart / 2)
    for k in gdat.indxmpolodim:
        indxmpol = np.where((gdat.meanmpol > gdat.binsmpolodim[k]) & (gdat.meanmpol < gdat.binsmpolodim[k+1]))
        psecodim[k] = np.mean(psec[indxmpol])
    psecodim *= gdat.meanmpolodim**2
    
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


def retr_eerrnorm(minmvarb, maxmvarb, meanvarb, stdvvarb):
   
    cdfnminm = 0.5 * (sp.special.erf((minmvarb - meanvarb) / stdvvarb / np.sqrt(2.)) + 1.)
    cdfnmaxm = 0.5 * (sp.special.erf((maxmvarb - meanvarb) / stdvvarb / np.sqrt(2.)) + 1.)
    cdfndiff = cdfnmaxm - cdfnminm
    
    return cdfnminm, cdfndiff
    

def retr_condcatl(gdat):
  
    if gdat.verbtype > 1:
        print 'listindxelemfull'
        print listindxelemfull
    
    # setup
    ## number of stacked samples
    numbstks = 0
    indxtupl = []
    indxstks = []
    indxstkssamp = []
    for n in gdat.indxsamptotl:
        indxstks.append([])
        indxstkssamptemp = []
        for l in gdat.fittindxpopl:
            indxstks[n].append([])
            for k in range(len(listindxelemfull[n][l])):
                indxstks[n][l].append(numbstks)
                indxstkssamptemp.append(numbstks)
                indxtupl.append([n, l, k])
                numbstks += 1
        indxstkssamp.append(np.array(indxstkssamptemp))
    
    if gdat.verbtype > 1:
        print 'indxstks'
        print indxstks
        print 'indxtupl'
        print indxtupl
        print 'indxstkssamp'
        print indxstkssamp
        print 'numbstks'
        print numbstks

    cntr = 0 
    arrystks = np.zeros((numbstks, gdat.fittnumbcomptotl))
    for n in gdat.indxsamptotl:
        indxsampcomp = retr_indxsampcomp(gdat, listindxelemfull[n], 'fitt') 
        for l in gdat.fittindxpopl:
            for k in np.arange(len(listindxelemfull[n][l])):
                for m, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                    arrystks[indxstks[n][l][k], m] = gdat.listsamp[n, indxsampcomp[strgcomp][l][k]]

    if gdat.verbtype > 0:
        print 'Constructing the distance matrix for %d stacked samples...' % arrystks.shape[0]
        timeinit = gdat.functime()
    
    gdat.distthrs = np.empty(gdat.fittnumbcomptotl)
    for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
        gdat.distthrs[k] = gdat.stdp[getattr(gdat, 'indxstdp' + strgcomp)]
    
    if gdat.verbtype > 0:
        print 'gdat.distthrs'
        print gdat.distthrs

    # construct lists of samples for each proposal type
    listdisttemp = [[] for k in range(gdat.fittnumbcomptotl)]
    indxstksrows = [[] for k in range(gdat.fittnumbcomptotl)]
    indxstkscols = [[] for k in range(gdat.fittnumbcomptotl)]
    thisperc = 0
    cntr = 0
    for k in gdat.fittindxcomptotl:
        for n in range(numbstks):
            dist = np.fabs(arrystks[n, k] - arrystks[:, k])
            indxstks = np.where(dist < gdat.distthrs[k])[0]
            if indxstks.size > 0:
                for j in indxstks:
                    cntr += 1
                    listdisttemp[k].append(dist[j])
                    indxstksrows[k].append(n)
                    indxstkscols[k].append(j)
            
            nextperc = np.floor(100. * float(k * numbstks + n) / numbstks / gdat.fittnumbcomptotl)
            if nextperc > thisperc:
                thisperc = nextperc
                print '%3d%% completed.' % thisperc
            if cntr > 1e6:
                break
        
        listdisttemp[k] = np.array(listdisttemp[k])
        indxstksrows[k] = np.array(indxstksrows[k])
        indxstkscols[k] = np.array(indxstkscols[k])

        if cntr > 1e6:
            break
    
    listdist = [[] for k in range(gdat.fittnumbcomptotl)]
    for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
        listdist[k] = scipy.sparse.csr_matrix((listdisttemp[k], (indxstksrows[k], indxstkscols[k])), shape=(numbstks, numbstks))
    
    listindxstkspair = []
    indxstksleft = []

    if gdat.verbtype > 0:
        timefinl = gdat.functime()
        print 'Done in %.3g seconds.' % (timefinl - timeinit)
    
    indxstksleft = range(numbstks)

    # list of sample lists of the labeled element
    indxstksassc = []
    cntr = 0
    
    gdat.prvlthrs = 0.05

    if gdat.verbtype > 1:
        print 'listdist'
        print listdist
    
    while len(indxstksleft) > 0:
        
        if gdat.verbtype > 1:
            print 'indxstksassc'
            print indxstksassc
            print 'indxstksleft'
            print indxstksleft
        
        # count number of associations
        numbdist = np.zeros(numbstks, dtype=int) - 1
        for p in range(len(indxstksleft)):
            indxindx = np.where((listdist[0][indxstksleft[p], :].tonp.array().flatten() * 2. * gdat.maxmlgal < gdat.anglassc) & \
                             (listdist[1][indxstksleft[p], :].tonp.array().flatten() * 2. * gdat.maxmbgal < gdat.anglassc))[0]
            numbdist[indxstksleft[p]] = indxindx.size
            
        prvlmaxmesti = np.amax(numbdist) / float(gdat.numbsamptotl)
        
        if prvlmaxmesti < gdat.prvlthrs:
            print 'Exiting...'
            break

        # determine the element with the highest number of neighbors
        indxstkscntr = np.argmax(numbdist)
        indxsamptotlcntr = indxtupl[indxstkscntr][0]
        indxpoplcntr = indxtupl[indxstkscntr][1]
        indxelemcntr = indxtupl[indxstkscntr][2]

        # add the central element sample
        indxstksassc.append([])
        indxstksassc[cntr].append(indxstkscntr)
        indxstksleft.remove(indxstkscntr)

        if gdat.verbtype > 1:
            print 'Match step %d' % cntr
            print 'numbdist'
            print numbdist
            print 'indxstkscntr'
            print indxstkscntr
            print 'indxstksleft'
            print indxstksleft
        
        # add the associated element samples
        if len(indxstksleft) > 0:
            for n in gdat.indxsamptotl:
                
                indxstkstemp = np.intersect1d(np.array(indxstksleft), indxstkssamp[n])
                
                if gdat.verbtype > 1:
                    print 'n'
                    print n
                    print 'indxstkstemp'
                    print indxstkstemp
                
                if n == indxsamptotlcntr:
                    continue
                
                if indxstkstemp.size > 0:
                    totl = np.zeros_like(indxstkstemp)
                    for k in gdat.fittindxcomptotl:
                        temp = listdist[k][indxstkscntr, indxstkstemp].tonp.array()[0]
                        totl = totl + temp**2

                    indxleft = np.argsort(totl)[0]
                    
                    if gdat.verbtype > 1:
                        print 'indxleft'
                        print indxleft
                    
                    indxstksthis = indxstkstemp[indxleft]
                
                    if gdat.verbtype > 1:
                        print 'indxstksthis'
                        print indxstksthis
                
                    thisbool = True
                    for k in gdat.fittindxcomptotl:
                        if listdist[k][indxstkscntr, indxstksthis] > gdat.distthrs[k]:
                            thisbool = False

                    if thisbool:
                        indxstksassc[cntr].append(indxstksthis)
                        indxstksleft.remove(indxstksthis)
                        if gdat.verbtype > 1:
                            print 'Appending...'
            
                # temp
                #if gdat.makeplot:
                #    gdatmodi = tdpy.util.gdatstrt()
                #    gdatmodi.thisindxelemfull = deepcopy(listindxelemfull[n])
                #    for r in range(len(indxstksassc)): 
                #        calc_poststkscond(gdat, indxstksassc)
                #    gdatmodi.thisindxelemfull = [[] for l in gdat.fittindxpopl]
                #    for indxstkstemp in indxstksleft:
                #        indxsamptotlcntr = indxtupl[indxstkstemp][0]
                #        indxpoplcntr = indxtupl[indxstkstemp][1]
                #        indxelemcntr = indxtupl[indxstkstemp][2]
                #        gdatmodi.thissamp = gdat.listsamp[indxsamptotlcntr, :]
                #        gdatmodi.thisindxelemfull[].append()

                #    plot_genemaps(gdat, gdatmodi, 'this', 'cntpdata', strgpdfn, indxenerplot=0, indxevttplot=0, cond=True)
                
                if gdat.verbtype > 1:
                    print 

            cntr += 1
        
        if gdat.verbtype > 1:
            print 
            print 
            print 
        
    if gdat.verbtype > 1:
        print 'listindxelemfull'
        print listindxelemfull
        print 'indxstksassc'
        print indxstksassc
        print
        for strgcomp in gdat.fittliststrgcomptotl:
            print 'strgcomp'
            print strgcomp
            print 'getattr(gdat, list + strgcomp)'
            print getattr(gdat, 'list' + strgcomp)

    gdat.dictglob['poststkscond'] = []
    gdat.dictglob['liststkscond'] = []
    # for each condensed element
    for r in range(len(indxstksassc)): 
        gdat.dictglob['liststkscond'].append([])
        gdat.dictglob['liststkscond'][r] = {}
        gdat.dictglob['poststkscond'].append([])
        gdat.dictglob['poststkscond'][r] = {}
        for strgfeat in gdat.fittliststrgfeattotl:
            gdat.dictglob['liststkscond'][r][strgfeat] = []

        # for each associated sample associated with the central stacked sample 
        for k in range(len(indxstksassc[r])):
            indxsamptotlcntr = indxtupl[indxstksassc[r][k]][0]
            indxpoplcntr = indxtupl[indxstksassc[r][k]][1]
            indxelemcntr = indxtupl[indxstksassc[r][k]][2]
            
            if gdat.verbtype > 1:
                print 'indxsamptotlcntr'
                print indxsamptotlcntr
                print 'indxpoplcntr'
                print indxpoplcntr
                print 'indxelemcntr'
                print indxelemcntr
                print
            
            for strgfeat in gdat.fittliststrgfeattotl:
                temp = getattr(gdat, 'list' + strgfeat)
                if temp[indxsamptotlcntr][indxpoplcntr].size > 0:
                    temp = temp[indxsamptotlcntr][indxpoplcntr][..., indxelemcntr]
                    gdat.dictglob['liststkscond'][r][strgfeat].append(temp)

    for r in range(len(gdat.dictglob['liststkscond'])):
        for strgfeat in gdat.fittliststrgfeattotl:
            arry = np.stack(gdat.dictglob['liststkscond'][r][strgfeat], axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat] = np.zeros(([3] + list(arry.shape[1:])))
            gdat.dictglob['poststkscond'][r][strgfeat][0, ...] = median(arry, axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat][1, ...] = percennp.tile(arry, 16., axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat][2, ...] = percennp.tile(arry, 84., axis=0)
            
    gdat.numbstkscond = len(gdat.dictglob['liststkscond'])

    if gdat.verbtype > 1:
        print 'gdat.dictglob[liststkscond]'
        for r in range(len(gdat.dictglob['liststkscond'])):
            print 'r'
            print r
            for strgfeat in gdat.fittliststrgfeatodimtotl:
                print strgfeat
                print gdat.dictglob['liststkscond'][r][strgfeat]
                print 
        print 'gdat.dictglob[poststkscond]'
        print gdat.dictglob['poststkscond']
    
    gdat.indxstkscond = np.arange(gdat.numbstkscond)
    gdat.prvl = np.empty(gdat.numbstkscond)
    for r in gdat.indxstkscond:
        gdat.prvl[r] = len(gdat.dictglob['liststkscond'][r]['deltllik'])
    gdat.prvl /= gdat.numbsamptotl
    gdat.minmprvl = 0.
    gdat.maxmprvl = 1.
    retr_axis(gdat, 'prvl', gdat.minmprvl, gdat.maxmprvl, 10)
    gdat.histprvl = np.histogram(gdat.prvl, bins=gdat.binsprvl)[0]
    if gdat.makeplot:
        pathcond = getattr(gdat, 'path' + strgpdfn + 'finlcond')
        for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
            path = pathcond + 'histdist' + strgcomp 
            listtemp = np.copy(listdist[k].tonp.array()).flatten()
            listtemp = listtemp[np.where(listtemp != 1e20)[0]]
            tdpy.mcmc.plot_hist(path, listtemp, r'$\Delta \tilde{' + getattr(gdat, 'labl' + strgcomp) + '}$')
            path = pathcond + 'histprvl'
            tdpy.mcmc.plot_hist(path, gdat.prvl, r'$p$')
    gdat.prvlthrs = 0.1 
    gdat.indxprvlhigh = np.where(gdat.prvl > gdat.prvlthrs)[0]
    gdat.numbprvlhigh = gdat.indxprvlhigh.size


def retr_conv(gdat, defl):
    
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    # temp
    conv = abs(np.gradient(defl[:, :, 0], gdat.sizepixl, axis=0) + np.gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) / 2.
    conv = conv.flatten()
    
    return conv


def retr_invm(gdat, defl):
    
    # temp
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    invm = (1. - np.gradient(defl[:, :, 0], gdat.sizepixl, axis=0)) * (1. - np.gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) - \
                                                np.gradient(defl[:, :, 0], gdat.sizepixl, axis=1) * np.gradient(defl[:, :, 1], gdat.sizepixl, axis=0)
    invm = invm.flatten()
    return invm


def setp_indxswepsave(gdat):

    gdat.indxswep = np.arange(gdat.numbswep)
    gdat.boolsave = np.zeros(gdat.numbswep, dtype=bool)
    gdat.indxswepsave = np.arange(gdat.numbburn, gdat.numbburn + gdat.numbsamp * gdat.factthin, gdat.factthin)
    gdat.boolsave[gdat.indxswepsave] = True
    gdat.indxsampsave = np.zeros(gdat.numbswep, dtype=int) - 1
    gdat.indxsampsave[gdat.indxswepsave] = np.arange(gdat.numbsamp)
    

def retr_cntspnts(gdat, listposi, spec):
    
    cnts = np.zeros((gdat.numbener, spec.shape[1]))
    
    if gdat.numbpixlfull > 1:
        lgal = listposi[0]
        bgal = listposi[1]
        indxpixlpnts = retr_indxpixl(gdat, bgal, lgal)
    else:
        elin = listposi[0]
        indxpixlpnts = np.zeros_like(elin, dtype=int)
    for k in range(spec.shape[1]):
        cnts[:, k] += spec[:, k] * gdat.expototl[:, indxpixlpnts[k]]
    if gdat.enerdiff:
        cnts *= gdat.deltener[:, None]
    cnts = np.sum(cnts, axis=0)

    return cnts


def retr_mdencrit(gdat, adissour, adishost, adishostsour):
    
    mdencrit = gdat.factnewtlght / 4. / np.pi * adissour / adishostsour / adishost
        
    return mdencrit


def retr_massfrombein(gdat, adissour, adishost, adishostsour):

    mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
    massfrombein = np.pi * adishost**2 * mdencrit

    return massfrombein


def retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut):
    
    mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
    
    fracacutasca = acut / asca
    
    factmcutfromdefs = np.pi * adishost**2 * mdencrit * asca * retr_mcutfrommscl(fracacutasca)

    return factmcutfromdefs


def retr_mcut(gdat, defs, asca, acut, adishost, mdencrit):
    
    mscl = defs * np.pi * adishost**2 * mdencrit * asca
    fracacutasca = acut / asca
    mcut = mscl * retr_mcutfrommscl(fracacutasca)
    
    return mcut


def retr_mcutfrommscl(fracacutasca):
    
    mcut = fracacutasca**2 / (fracacutasca**2 + 1.)**2 * ((fracacutasca**2 - 1.) * np.log(fracacutasca) + fracacutasca * np.pi - (fracacutasca**2 + 1.))

    return mcut


def retr_negalogt(varb):
    
    negalogt = sign(varb) * np.log10(np.fabs(varb))
    
    return negalogt


def setp_prem(gdat):

    # set up the indices of the models
    for strgmodl in gdat.liststrgmodl:
        retr_indxsamp(gdat, strgmodl=strgmodl, init=True)
    
    gdat.liststrglimt = ['minm', 'maxm']
   
    if gdat.numbswepplot is None:
        gdat.numbswepplot = 50000
   
    gdat.numbplotfram = gdat.numbswep / gdat.numbswepplot

    gdat.refrcolr = 'mediumseagreen'
    gdat.fittcolr = 'deepskyblue'

    gdat.minmmass = 1.
    gdat.maxmmass = 10.
    
    if gdat.checprio:
        gdat.liststrgpdfn = ['prio', 'post']
    else:
        gdat.liststrgpdfn = ['post']

    gdat.lablmass = 'M'
    gdat.minmmassshel = 1e1
    gdat.maxmmassshel = 1e5
    gdat.lablmassshel = '$M_r$' 

    gdat.lablcurv = r'\kappa'
    gdat.lablexpc = r'E_{c}'
    
    gdat.scalcurvplot = 'self'
    gdat.scalexpcplot = 'self'
    gdat.factcurvplot = 1.
    gdat.factexpcplot = 1.
    
    gdat.minmper0 = 1e-3 
    gdat.maxmper0 = 1e1
    
    gdat.minmmagf = 10**7.5
    gdat.maxmmagf = 10**16
    
    # temp -- automatize this eventually
    gdat.trueminmper0 = gdat.minmper0
    gdat.fittminmper0 = gdat.minmper0
    gdat.truemaxmper0 = gdat.maxmper0
    gdat.fittmaxmper0 = gdat.maxmper0
    gdat.trueminmmagf = gdat.minmmagf
    gdat.fittminmmagf = gdat.minmmagf
    gdat.truemaxmmagf = gdat.maxmmagf
    gdat.fittmaxmmagf = gdat.maxmmagf

    gdat.fittlistelemmrkr = ['+', '_', '3']
    gdat.refrlistelemmrkrhits = ['x', '|', '4']
    gdat.refrlistelemmrkrmiss = ['s', 'o', 'p']
    
    gdat.listscaltype = ['self', 'logt', 'atan', 'gaus', 'pois']
    
    gdat.numbgrid = 1
    gdat.indxgrid = np.arange(gdat.numbgrid)

    if gdat.pixltype == 'heal' and gdat.forccart:
        raise Exception('Cartesian forcing can only used with cart pixltype')

    gdat.liststrgphas = ['fram', 'finl', 'anim']
    gdat.liststrgelemtdimtype = ['bind']
    
    # input data
    if gdat.datatype == 'inpt':
        path = gdat.pathinpt + gdat.strgexprsbrt
        gdat.sbrtdata = astropy.io.fits.getdata(path)
            
        if gdat.pixltype == 'heal' or gdat.pixltype == 'cart' and gdat.forccart:
            if gdat.sbrtdata.ndim != 3:
                raise Exception('exprsbrtdata should be a 3D numpy np.array if pixelization is HealPix.')
        else:
            if gdat.sbrtdata.ndim != 4:
                print 'gdat.sbrtdata'
                summgene(gdat.sbrtdata)
                raise Exception('exprsbrtdata should be a 4D numpy np.array if pixelization is Cartesian.')
        
        if gdat.pixltype == 'cart' and not gdat.forccart:
            gdat.sbrtdata = gdat.sbrtdata.reshape((gdat.sbrtdata.shape[0], -1, gdat.sbrtdata.shape[3]))
                    
        gdat.numbenerfull = gdat.sbrtdata.shape[0]
        if gdat.pixltype == 'heal':
            gdat.numbpixlfull = gdat.sbrtdata.shape[1]
        elif gdat.forccart:
            gdat.numbpixlfull = gdat.numbsidecart**2
        else:
            gdat.numbpixlfull = gdat.sbrtdata.shape[1] * gdat.sbrtdata.shape[2]
        gdat.numbevttfull = gdat.sbrtdata.shape[2]
        
        if gdat.pixltype == 'heal':
            # temp
            gdat.numbsidecart = 100
            gdat.numbsideheal = int(np.sqrt(gdat.numbpixlfull / 12))
    

def setpinit(gdat, boolinitsetp=False):
    
    if False and gdat.boolelemdeflsubhanyy and gdat.strgproc == 'fink2.rc.fas.harvard.edu':
        cliblens = ctypes.CDLL(os.environ["PCAT_PATH"] + '/cliblens.so')
        cliblens.retr_deflsubh()
    
    # set up the indices of the fitting model
    retr_indxsamp(gdat)
    
    gdat.refrliststrgfeattagg = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
    for q in gdat.indxrefr:
        for strgfeat in gdat.refrliststrgfeat[q]:
            for l in gdat.fittindxpopl:
                if strgfeat in gdat.refrliststrgfeatonly[q][l]:
                    gdat.refrliststrgfeattagg[q][l].append(strgfeat + gdat.listnamerefr[q])
                else:
                    gdat.refrliststrgfeattagg[q][l].append(strgfeat)
    
    # common element features
    gdat.liststrgfeatcomm = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
    for q in gdat.indxrefr:
        for l in gdat.fittindxpopl:
            for strgfeat in gdat.fittliststrgfeat[l]:
                if strgfeat in gdat.refrliststrgfeat[q]:
                    gdat.liststrgfeatcomm[q][l].append(strgfeat)
    
    for l in gdat.fittindxpopl:
        for strgfeat in gdat.fittliststrgfeat[l]:
            if strgfeat.startswith('aerr'):
                setattr(gdat, 'minm' + strgfeat, -100.)
                setattr(gdat, 'maxm' + strgfeat, 100.)

    for strgvarb in ['boolelempsfnanyy']:
        varbcomm = False
        for strgmodl in gdat.liststrgmodl:
            varb = getattr(gdat, strgmodl + strgvarb)
            varbcomm = varbcomm or varb
        setattr(gdat, 'comm' + strgvarb, varbcomm) 

    gdat.listnamevarbstat = ['samp', 'sampunit', 'indxelemfull', 'lliktotl', 'llik', 'lpritotl', 'lpri']
    if gdat.pixltype == 'cart' and (gdat.fittpsfnevaltype == 'conv' or gdat.fittpsfnevaltype == 'full'):
        gdat.listnamevarbstat += ['psfnconv']
    if gdat.fittboolelemsbrtdfncanyy:
        gdat.listnamevarbstat += ['sbrtdfnc']
    if gdat.fittboolelemsbrtextsbgrdanyy:
        gdat.listnamevarbstat += ['sbrtextsbgrd']
    if gdat.fittlensmodltype != 'none':
        gdat.listnamevarbstat += ['sbrtlens']
    if gdat.fittlensmodltype != 'none' or gdat.fitthostemistype != 'none':
        for e in gdat.fittindxsersfgrd:
            if gdat.fittlensmodltype != 'none':
                gdat.listnamevarbstat += ['deflhostisf%d' % e]
            if gdat.fitthostemistype != 'none':
                gdat.listnamevarbstat += ['sbrthostisf%d' % e]
    if gdat.fittconvdiffanyy and (gdat.fittpsfnevaltype == 'full' or gdat.fittpsfnevaltype == 'conv'):
        gdat.listnamevarbstat += ['sbrtmodlconv']
    if gdat.fittboolelemdeflsubhanyy:
        gdat.listnamevarbstat += ['deflsubh']
    
    # paths
    ## data
    gdat.pathpixlcnvt = gdat.pathdata + 'pixlcnvt/'
    gdat.pathprox = gdat.pathdata + 'prox/'
    ## plot
    gdat.pathplotrtag = gdat.pathimag + gdat.rtag + '/'
    gdat.pathinit = gdat.pathplotrtag + 'init/'
    gdat.pathinitintr = gdat.pathinit + 'intr/'
    
    ## names of the variables for which cumulative posteriors will be plotted
    if 'lens' in gdat.fittelemtype:
        gdat.listnamevarbcpct = ['convelem']
    else:
        gdat.listnamevarbcpct = []
    
    if gdat.numbpixlfull != 1:
        gdat.ascaglob = 0.05 / gdat.anglfact
        gdat.acutglob = 1. / gdat.anglfact
        gdat.cutfdefs = 3e-3 / gdat.anglfact

    # plotting
    gdat.legdsampdist = 'Posterior'
    gdat.legdsamp = 'Sample'
    gdat.legdmlik = 'Maximum likelihood'
    gdat.legdmedi = 'Median'
    gdat.legdpmea = 'Mean'
    gdat.legdstdv = 'Std. dev.'
  
    # number of samples for which cumulative posterior will be calculated
    gdat.numbsampcpct = 10
    gdat.indxsampcpct = np.arange(gdat.numbsampcpct)
    
    # p value contours 
    gdat.pvalcont = [0.317, 0.0455, 2.7e-3, 6e-5, 1.3e-6]

    ## number of bins in histogram plots
    gdat.numbbinsplot = 20
    gdat.indxbinsplot = np.arange(gdat.numbbinsplot)
    
    ## number of bins in hyperprior plots
    gdat.numbbinsplotprio = 100
    # temp
    ### number of standard deviations away from the mean of spectral parameters to plot
    gdat.numbstdvspepdist = 3.
    ### minima and maxima for spectral parameters
    gdat.numbstdv = 3.
    
    ## labels and scales for variables
    if gdat.fittlensmodltype != 'none':
        setattr(gdat, 'lablmasssubhintg', r'$M_{\rm{sub}}$')
        setattr(gdat, 'lablmasssubhdelt', r'$\rho_{\rm{sub}}$')
        setattr(gdat, 'lablmasssubhintgbein', r'$M_{\rm{sub,E}}$')
        setattr(gdat, 'lablmasssubhdeltbein', r'$\rho_{\rm{sub,E}}$')
        setattr(gdat, 'lablmasssubhintgunit', '$10^9 M_{\odot}$')
        setattr(gdat, 'lablmasssubhdeltunit', '$M_{\odot}$/kpc')
        setattr(gdat, 'factmasssubhintgplot', 1e-9)
        setattr(gdat, 'factmasssubhintgbeinplot', 1e-9)
        setattr(gdat, 'lablmasssubhintgbeinunit', '$10^9 M_{\odot}$')
        setattr(gdat, 'lablmasssubhdeltbeinunit', '$M_{\odot}$/kpc')
        setattr(gdat, 'lablfracsubhintg', r'f_{\rm{sub}}')
        setattr(gdat, 'lablfracsubhdelt', r'f_{\rho,\rm{sub}}')
        setattr(gdat, 'lablfracsubhintgbein', r'$f_{\rm{sub,E}}$')
        setattr(gdat, 'lablfracsubhdeltbein', r'$f_{\rho,\rm{sub,E}}$')
        for e in gdat.fittindxsersfgrd:
            setattr(gdat, 'lablmasshostisf%dbein' % e, r'$M_{\rm{hst,%d,C}}$' % e)
            setattr(gdat, 'lablmasshostisf%dintg' % e, r'$M_{\rm{hst,%d<}}$' % e)
            setattr(gdat, 'lablmasshostisf%ddelt' % e, r'$M_{\rm{hst,%d}}$' % e)
            setattr(gdat, 'lablmasshostisf%dintgbein' % e, r'$M_{\rm{hst,E,%d<}}$' % e)
            setattr(gdat, 'lablmasshostisf%ddeltbein' % e, r'$M_{\rm{hst,E,%d}}$' % e)
        for namevarb in ['fracsubh', 'masssubh']:
            for namecalc in ['delt', 'intg']:
                for nameeval in ['', 'bein']:
                    setattr(gdat, 'scal' + namevarb + namecalc + nameeval, 'logt')
        for e in gdat.fittindxsersfgrd:
            setattr(gdat, 'scalmasshostisf%d' % e + 'bein', 'logt')
            for namecalc in ['delt', 'intg']:
                for nameeval in ['', 'bein']:
                    setattr(gdat, 'scalmasshostisf%d' % e + namecalc + nameeval, 'logt')
    
    # scalar variable setup
    gdat.lablhistcntplowrdfncsubten00evt0 = 'N_{pix,l}'
    gdat.lablhistcntphigrdfncsubten00evt0 = 'N_{pix,h}'
    gdat.lablhistcntplowrdfncen00evt0 = 'N_{pix,l}'
    gdat.lablhistcntphigrdfncen00evt0 = 'N_{pix,h}'
    for i in gdat.indxener:
        setattr(gdat, 'lablfracsdenmeandarkdfncsubten%02d' % i, 'f_{D/ST,%d}' % i)
    gdat.lablbooldfncsubt = 'H'
    
    gdat.lablpriofactdoff = r'$\alpha_{p}$'
    gdat.scalpriofactdoff = 'self'

    gdat.minmreds = 0.
    gdat.maxmreds = 1.5
    
    gdat.minmmagt = 19.
    gdat.maxmmagt = 28.

    gdat.scalmaxmnumbelem = 'logt'
    gdat.scalmedilliktotl = 'logt'

    gdat.lablener = 'E'
    #gdat.lablenertotl = '$%s$ [%s]' % (gdat.lablener, gdat.strgenerunit)
    
    gdat.lablgwdt = r'\sigma_G'
    
    gdat.lablgang = r'\theta'
    gdat.lablaang = r'\phi'
    gdat.labllgalunit = gdat.lablgangunit
    gdat.lablbgalunit = gdat.lablgangunit
   
    gdat.lablanglfromhost = r'\theta_{\rm{0,hst}}'
    gdat.lablanglfromhostunit = gdat.lablgangunit

    gdat.labldefs = r'\alpha_s'
    gdat.lablflux = 'f'
    gdat.lablnobj = 'p'
    
    gdat.lablelin = r'\mathcal{E}'
    
    gdat.lablsbrt = r'\Sigma'
    
    gdat.labldeflprof = r'\alpha_a'
    gdat.labldeflprofunit = u'$^{\prime\prime}$'
    
    gdat.strgenerkevv = 'keV'
    gdat.strgenergevv = 'GeV'
    gdat.strgenerergs = 'erg'
    gdat.strgenerimum = '\mu m^{-1}'

    gdat.labldefsunit = u'$^{\prime\prime}$'
    gdat.lablprat = 'cm$^{-2}$ s$^{-1}$'
    for nameenerscaltype in ['en00', 'en01', 'en02', 'en03']:
        
        for labltemptemp in ['flux', 'sbrt']:
            labltemp = getattr(gdat, 'labl' + labltemptemp)

            # define the label
            if nameenerscaltype == 'en00':
                strgenerscal = '%s' % labltemp
            if nameenerscaltype == 'en01':
                strgenerscal = 'E%s' % labltemp
            if nameenerscaltype == 'en02':
                strgenerscal = 'E^2%s' % labltemp
            if nameenerscaltype == 'en03':
                strgenerscal = '%s' % labltemp
            labl = '%s' % strgenerscal
            setattr(gdat, 'labl' + labltemptemp + nameenerscaltype, labl)

            for nameenerunit in ['gevv', 'ergs', 'kevv', 'imum']:
                
                strgenerunit = getattr(gdat, 'strgener' + nameenerunit)

                if nameenerscaltype == 'en00':
                    strgenerscalunit = '%s$^{-1}$' % strgenerunit
                if nameenerscaltype == 'en01':
                    strgenerscalunit = '' 
                if nameenerscaltype == 'en02':
                    strgenerscalunit = '%s' % strgenerunit
                if nameenerscaltype == 'en03':
                    strgenerscalunit = '%s' % strgenerunit
                
                # define the label unit
                for namesoldunit in ['ster', 'degr']:
                    if labltemptemp == 'flux':
                        lablunit = '%s %s' % (strgenerscalunit, gdat.lablprat)
                        setattr(gdat, 'lablflux' + nameenerscaltype + nameenerunit + 'unit', lablunit)
                    else:
                        if namesoldunit == 'ster':
                            lablunit = '%s %s sr$^{-1}$' % (strgenerscalunit, gdat.lablprat)
                        if namesoldunit == 'degr':
                            lablunit = '%s %s deg$^{-2}$' % (strgenerscalunit, gdat.lablprat)
                        setattr(gdat, 'lablsbrt' + nameenerscaltype + nameenerunit + namesoldunit + 'unit', lablunit)

                    # define the total label
                    if labltemptemp == 'flux':
                        setattr(gdat, 'lablflux' + nameenerscaltype + nameenerunit + 'totl', '$%s$ [%s]' % (labl, lablunit))
                    else:
                        setattr(gdat, 'lablsbrt' + nameenerscaltype + nameenerunit + namesoldunit + 'totl', '$%s$ [%s]' % (labl, lablunit))
    if gdat.enerbins:
        gdat.lablfluxunit = getattr(gdat, 'lablfluxen00' + gdat.nameenerunit + 'unit')
        gdat.lablsbrtunit = getattr(gdat, 'lablsbrten00' + gdat.nameenerunit + 'sterunit')

    gdat.lablexpo = r'$\epsilon$'
    gdat.lablexpounit = 'cm$^2$ s'
    
    gdat.lablprvl = '$p$'
    
    gdat.lablreds = 'z'
    gdat.lablmagt = 'm_R'
    
    setattr(gdat, 'lablper0', 'P_0')
    setattr(gdat, 'factper0plot', 1.)
    setattr(gdat, 'scalper0plot', 'logt')
  
    gdat.labldglc = 'd_{gc}'
    gdat.factdglcplot = 1e-3
    gdat.scaldglcplot = 'logt'
    
    gdat.labldlos = 'd_{los}'
    gdat.scaldlosplot = 'logt'
    if gdat.exprtype == 'ferm':
        gdat.labldlosunit = 'kpc'
        gdat.labllumi = r'L_{\gamma}'
    if gdat.exprtype == 'chan':
        gdat.labldlosunit = 'Mpc'
        gdat.labllumi = r'L_{X}'
        gdat.labllum0 = r'L_{X, 0}'
    
    gdat.lablgeff = r'\eta_{\gamma}'
    gdat.factgeffplot = 1.
    gdat.scalgeffplot = 'logt'
    
    gdat.factlumiplot = 1.
    gdat.scallumiplot = 'logt'
    gdat.labllumiunit = 'erg s$^{-1}$'
    gdat.labllum0unit = 'erg s$^{-1}$'
    
    gdat.lablthet = r'\theta_{gc}'
    gdat.factthetplot = 180. / np.pi
    gdat.scalthetplot = 'self'
    
    gdat.lablphii = r'\phi_{gc}'
    gdat.factphiiplot = 1.
    gdat.scalphiiplot = 'self'
    
    setattr(gdat, 'lablmagf', 'B')
    setattr(gdat, 'factmagfplot', 1.)
    setattr(gdat, 'scalmagfplot', 'logt')
    
    setattr(gdat, 'lablper1', 'P_1')
    if gdat.datatype == 'inpt':
        setattr(gdat, 'minmper0', 1e-3)
        setattr(gdat, 'maxmper0', 1e1)
        setattr(gdat, 'minmper1', 1e-20)
        setattr(gdat, 'maxmper1', 1e-10)
        setattr(gdat, 'minmper1', 1e-20)
        setattr(gdat, 'maxmper1', 1e-10)
        setattr(gdat, 'minmflux0400', 1e-1)
        setattr(gdat, 'maxmflux0400', 1e4)
    setattr(gdat, 'factper1plot', 1.)
    setattr(gdat, 'scalper1plot', 'logt')
    setattr(gdat, 'lablflux0400', 'S_{400}')
    setattr(gdat, 'factflux0400plot', 1.)
    setattr(gdat, 'scalflux0400plot', 'logt')
    
    for q in gdat.indxrefr:
        setattr(gdat, 'lablaerr' + gdat.listnamerefr[q], '\Delta_{%d}' % q)
    gdat.lablsigm = '\sigma_l'
    gdat.lablgamm = '\gamma_l'

    gdat.lablbcom = '\eta'
    
    gdat.lablinfopost = 'D_{KL}'
    gdat.lablinfopostunit = 'nat'
    gdat.lablinfoprio = 'D_{KL,pr}'
    gdat.lablinfopriounit = 'nat'
    
    gdat.labllevipost = '\ln P(D)'
    gdat.labllevipostunit = 'nat'
    gdat.lablleviprio = '\ln P_{pr}(D)'
    gdat.labllevipriounit = 'nat'
    
    gdat.lablsind = 's'
    for i in gdat.indxenerinde:
        setattr(gdat, 'lablsindcolr%04d' % i, 's_%d' % i)

    gdat.lablexpcunit = gdat.strgenerunit
    
    gdat.labllliktotl = r'\ln P(D|M)'
    
    gdat.labllpripena = r'\ln P(N)'
    
    gdat.lablasca = r'\theta_s'
    gdat.lablascaunit = gdat.lablgangunit
    gdat.lablacut = r'\theta_c'
    gdat.lablacutunit = gdat.lablgangunit
    
    gdat.lablmcut = r'M_{c,n}'
    gdat.lablmcutunit = r'$M_{\odot}$'
    
    gdat.lablmcutcorr = r'\bar{M}_{c,n}'
    gdat.lablmcutcorrunit = r'$M_{\odot}$'
    
    gdat.lablspec = gdat.lablflux
    gdat.lablspecunit = gdat.lablfluxunit
    gdat.lablspecplot = gdat.lablflux
    gdat.lablspecplotunit = gdat.lablfluxunit
    gdat.lablcnts = 'C'
    gdat.labldeltllik = r'\Delta_n \ln P(D|M)'
    gdat.labldiss = r'\theta_{sa}'
    gdat.labldissunit = gdat.lablgangunit
    
    gdat.lablrele = r'\langle|\vec{\alpha}_n \cdot \vec{\nabla} k_l| \rangle'
    
    gdat.lablrelc = r'\langle\vec{\alpha}_n \cdot \vec{\nabla} k_l \rangle'
    
    gdat.lablreld = r'\langle|\vec{\alpha}_n \cdot \vec{\nabla} k_d| \rangle'
    
    gdat.lablreln = r'\langle \Delta \theta_{pix} |\hat{\alpha}_n \cdot \vec{\nabla} k_l| / \alpha_{s,n} \rangle'
    
    gdat.lablrelm = r'\langle |\vec{\nabla}_{\hat{\alpha}} k_l| / \alpha_{s,n} \rangle'
    gdat.lablrelk = r'\langle |\vec{\nabla}_{\hat{\alpha}} k_l| / \alpha_{s,n} \rangle'
    gdat.lablrelf = r'\langle |\vec{\nabla}_{\hat{\alpha}} k_l| / \alpha_{s,n} \rangle / k_m'
    
    gdat.minmcmpl = -0.2
    gdat.maxmcmpl = 1.2
    gdat.minmfdis = -0.2
    gdat.maxmfdis = 1.2
   
    for q in gdat.indxrefr:
        for l in gdat.fittindxpopl:
            setattr(gdat, 'minmcmplpop%dpop%d' % (l, q), gdat.minmcmpl)
            setattr(gdat, 'maxmcmplpop%dpop%d' % (l, q), gdat.maxmcmpl)
            setattr(gdat, 'corrcmplpop%dpop%d' % (l, q), None)
            setattr(gdat, 'factcmplpop%dpop%dplot' % (l, q), 1.)
            setattr(gdat, 'scalcmplpop%dpop%d' % (l, q), 'self')
            setattr(gdat, 'lablcmplpop%dpop%d' % (l, q), '$c_{%d%d}$' % (l, q))
            
            setattr(gdat, 'minmfdispop%dpop%d' % (q, l), gdat.minmfdis)
            setattr(gdat, 'maxmfdispop%dpop%d' % (q, l), gdat.maxmfdis)
            setattr(gdat, 'corrfdispop%dpop%d' % (q, l), None)
            setattr(gdat, 'factfdispop%dpop%dplot' % (q, l), 1.)
            setattr(gdat, 'scalfdispop%dpop%d' % (q, l), 'self')
            setattr(gdat, 'lablfdispop%dpop%d' % (q, l), '$f_{%d%d}$' % (q, l))
                
    dicttemp = deepcopy(gdat.__dict__)
    for name, valu in dicttemp.iteritems():
        if name.startswith('labl') and not name.endswith('unit'):
            name = name[4:]
            labl = getattr(gdat, 'labl' + name)
            try:
                lablunit = getattr(gdat, 'labl' + name + 'unit')
                if lablunit == '' or lablunit is None:
                    raise
                if labl.startswith('$'):
                    setattr(gdat, 'labl' + name + 'totl', '%s [%s]' % (labl, lablunit))
                else:
                    setattr(gdat, 'labl' + name + 'totl', '$%s$ [%s]' % (labl, lablunit))
            except:
                lablunit = ''
                setattr(gdat, 'labl' + name + 'unit', lablunit)
                if labl.startswith('$'):
                    setattr(gdat, 'labl' + name + 'totl', '%s' % labl)
                else:
                    setattr(gdat, 'labl' + name + 'totl', '$%s$' % labl)
    
            for strgextn in ['ref', 'sam']:
                if strgextn == 'ref':
                    nameextn = 'refr'
                if strgextn == 'sam':
                    nameextn = 'samp'
                lablxdat = '$%s^{%s}$%s' % (labl, strgextn, lablunit)
                if lablunit != '':
                    setattr(gdat, 'labl' + name + nameextn, '$%s^{%s}$ [%s]' % (labl, strgextn, lablunit))
                else:
                    setattr(gdat, 'labl' + name + nameextn, '$%s^{%s}$' % (labl, strgextn))
    
    if gdat.exprtype == 'chan':
        if gdat.anlytype == 'spec':
            gdat.minmspec = 1e-2
            gdat.maxmspec = 1e1
        else:
            gdat.minmspec = 1e-11
            gdat.maxmspec = 1e-7
    else:
        gdat.minmspec = 1e-11
        gdat.maxmspec = 1e-7
    
    if gdat.exprtype == 'ferm':
        gdat.minmlumi = 1e32
        gdat.maxmlumi = 1e36
    elif gdat.exprtype == 'chan':
        if gdat.datatype == 'inpt':
            gdat.minmlum0 = 1e42
            gdat.maxmlum0 = 1e46
        gdat.minmlumi = 1e41
        gdat.maxmlumi = 1e45
    
    try:
        gdat.minmdlos
    except:
        if gdat.exprtype == 'chan':
            gdat.minmdlos = 1e7
            gdat.maxmdlos = 1e9
        else:
            gdat.minmdlos = 6e3
            gdat.maxmdlos = 1.1e4
    
    if gdat.exprtype == 'ferm':
        gdat.minmcnts = 1e1
        gdat.maxmcnts = 1e5
    if gdat.exprtype == 'chan':
        if gdat.numbpixlfull == 1:
            gdat.minmcnts = 1e4
            gdat.maxmcnts = 1e8
        else:
            gdat.minmcnts = 1.
            gdat.maxmcnts = 1e3
    if gdat.exprtype == 'hubb':
        gdat.minmcnts = 1.
        gdat.maxmcnts = 1e3
    if gdat.exprtype == 'fire':
        gdat.minmcnts = 1.
        gdat.maxmcnts = 1e3

    gdat.minmspecplot = gdat.minmspec
    gdat.maxmspecplot = gdat.maxmspec
    
    gdat.minmdeltllik = 1.
    gdat.maxmdeltllik = 1e3
    gdat.minmdiss = 0.
    gdat.maxmdiss = gdat.maxmgangdata * np.sqrt(2.)
    
    gdat.minmrele = 1e-3
    gdat.maxmrele = 1e1

    gdat.minmreln = 1e-3
    gdat.maxmreln = 1.

    gdat.minmrelk = 1e-3
    gdat.maxmrelk = 1.

    gdat.minmrelf = 1e-5
    gdat.maxmrelf = 1e-1

    gdat.minmrelm = 1e-3
    gdat.maxmrelm = 1e1

    gdat.minmreld = 1e-3
    gdat.maxmreld = 1e1

    gdat.minmrelc = 1e-3
    gdat.maxmrelc = 1.

    gdat.minmmcut = 3e7
    gdat.maxmmcut = 2e9
    gdat.minmmcutcorr = gdat.minmmcut
    gdat.maxmmcutcorr = gdat.maxmmcut

    if gdat.numbpixlfull != 1:
        gdat.minmbein = 0.
        gdat.maxmbein = 1. / gdat.anglfact
    
    # scalar variables
    if gdat.numbpixlfull != 1:
        gdat.minmdeflprof = 1e-3 / gdat.anglfact
        gdat.maxmdeflprof = 0.1 / gdat.anglfact
    
    gdat.minmfracsubh = 0.
    gdat.maxmfracsubh = 0.3
    gdat.scalfracsubh = 'self'

    gdat.minmmasshost = 1e10
    gdat.maxmmasshost = 1e13
    gdat.scalmasshost = 'self'
    
    gdat.minmmasssubh = 1e8
    gdat.maxmmasssubh = 1e10
    gdat.scalmasssubh = 'self'

    if gdat.datatype == 'inpt':
        for l in gdat.fittindxpopl:
            for strgpdfn in gdat.fittliststrgpdfnprio[l]:
                if strgpdfn.startswith('gaum') and gdat.fittlgalprio is None and gdat.fittbgalprio is None:
                    raise Exception('If spatdisttype is "gaus", spatial coordinates of the prior catalog should be provided via lgalprio and bgalprio.')
    
    # temp -- have these definitions separate for all features
    # feature plotting factors and scalings
    gdat.lablfeat = {}
    gdat.dictglob = {}
    gdat.lablfeatunit = {}
    gdat.lablfeattotl = {}
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        for l in indxpopl:
            for strgfeat in liststrgfeat[l]:
                if strgfeat.startswith('defs') or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal' or strgfeat == 'gwdt' or \
                                                                                       strgfeat == 'diss' or strgfeat == 'asca' or strgfeat == 'acut':
                    setattr(gdat, 'fact' + strgfeat + 'plot', gdat.anglfact)
                else:
                    setattr(gdat, 'fact' + strgfeat + 'plot', 1.)
                if strgfeat == 'flux' or strgfeat == 'expo' or strgfeat == 'magf' or strgfeat == 'nobj' or \
                    strgfeat == 'relk' or strgfeat == 'relf' or strgfeat == 'elin' or strgfeat == 'flux0400' or \
                    strgfeat == 'cnts' or strgfeat.startswith('per') or strgfeat == 'gwdt' or strgfeat == 'dglc' or strgfeat.startswith('lumi') or \
                    strgfeat == 'relm' or strgfeat.startswith('dlos') or strgfeat.startswith('lum0') or strgfeat == 'defs' or \
                    strgfeat == 'relc' or strgfeat == 'rele' or strgfeat == 'reln' or strgfeat == 'reld' \
                    or strgfeat.startswith('mcut') or strgfeat == 'deltllik':
                    setattr(gdat, 'scal' + strgfeat + 'plot', 'logt')
                else:
                    setattr(gdat, 'scal' + strgfeat + 'plot', 'self')
    
    if gdat.exprtype == 'ferm':
        gdat.factdlosplot = 1e-3
    if gdat.exprtype == 'chan':
        gdat.factdlosplot = 1e-6

    # construct the fitting model
    setp_fixp(gdat)
    
    if gdat.datatype == 'mock':
        setp_fixp(gdat, strgmodl='true')
    
    for strgmodl in gdat.liststrgmodl:
        namefixp = getattr(gdat, strgmodl + 'namefixp')
        lablfixp = getattr(gdat, strgmodl + 'lablfixp')
        lablfixpunit = getattr(gdat, strgmodl + 'lablfixpunit')
        numbfixp = getattr(gdat, strgmodl + 'numbfixp')
        scalfixp = getattr(gdat, strgmodl + 'scalfixp')
        lablfixptotl = np.empty(numbfixp, dtype=object)
        for k in getattr(gdat, strgmodl + 'indxfixp'):
            if lablfixpunit[k] == '':
                lablfixptotl[k] = '%s' % lablfixp[k]
            else:
                lablfixptotl[k] = '%s [%s]' % (lablfixp[k], lablfixpunit[k])
        setattr(gdat, strgmodl + 'lablfixptotl', lablfixptotl)
        for k, name in enumerate(namefixp):
            setattr(gdat, 'labl' + namefixp[k] + 'totl', lablfixptotl[k])
            setattr(gdat, 'scal' + namefixp[k], scalfixp[k])
    
    if gdat.datatype != 'mock':
        if gdat.fittboolelemlghtanyy and gdat.exprtype == 'ferm' and gdat.maxmgangdata == 20. / gdat.anglfact:
            path = gdat.pathinpt + 'sbrt0018.png'
            gdat.sbrt0018 = sp.ndimage.imread(path, flatten=True)
            gdat.sbrt0018 -= np.amin(gdat.sbrt0018)
            gdat.sbrt0018 /= np.amax(gdat.sbrt0018)
            binslgaltemp = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.sbrt0018.shape[1])
            binsbgaltemp = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.sbrt0018.shape[0])
            gdat.sbrt0018objt = sp.interpolate.RectBivariateSpline(binsbgaltemp, binslgaltemp, gdat.sbrt0018)

    # list of scalar variable names
    gdat.listnamevarbscal = list(gdat.fittnamefixp) 
    booltemp = False
    if gdat.fittlensmodltype != 'none':
        for e in gdat.fittindxsersfgrd:
            strgisfr = 'isf%d' % e
            gdat.listnamevarbscal += ['masshost' + strgisfr + 'bein']
            setattr(gdat, 'minmmasshost' + strgisfr + 'bein', gdat.minmmasshost)
            setattr(gdat, 'maxmmasshost' + strgisfr + 'bein', gdat.maxmmasshost)
            for strgtemp in ['delt', 'intg']:
                gdat.listnamevarbscal += ['masshost' + strgisfr + strgtemp + 'bein']
                setattr(gdat, 'minmmasshost' + strgisfr + strgtemp + 'bein', gdat.minmmasshost)
                setattr(gdat, 'maxmmasshost' + strgisfr + strgtemp + 'bein', gdat.maxmmasshost)
        if gdat.fittnumbtrap > 0:
            if 'lens' in gdat.fittelemtype:
                for strgtemp in ['delt', 'intg']:
                    gdat.listnamevarbscal += ['masssubh' + strgtemp + 'bein', 'fracsubh' + strgtemp + 'bein'] 
                    setattr(gdat, 'minmmasssubh' + strgtemp + 'bein', gdat.minmmasssubh)
                    setattr(gdat, 'maxmmasssubh' + strgtemp + 'bein', gdat.maxmmasssubh)
                    setattr(gdat, 'minmfracsubh' + strgtemp + 'bein', gdat.minmfracsubh)
                    setattr(gdat, 'maxmfracsubh' + strgtemp + 'bein', gdat.maxmfracsubh)
    gdat.listnamevarbscal += ['lliktotl']
    if gdat.fittnumbtrap > 0:
        gdat.listnamevarbscal += ['lpripena']
    
    if gdat.fittboolelemsbrtdfncanyy:
        for strgbins in ['lowr', 'higr']:
            gdat.listnamevarbscal += ['histcntp%sdfncen00evt0' % strgbins]
            gdat.listnamevarbscal += ['histcntp%sdfncsubten00evt0' % strgbins]
        for i in gdat.indxener:
            gdat.listnamevarbscal += ['fracsdenmeandarkdfncsubten%02d' % i]
        gdat.listnamevarbscal += ['booldfncsubt']
            
    gdat.numbvarbscal = len(gdat.listnamevarbscal)
    gdat.indxvarbscal = np.arange(gdat.numbvarbscal)
        
    # plotting factors for scalar variables
    for name in gdat.listnamevarbscal:
        if name in gdat.fittnamefixp:
            indxfixptemp = np.where(gdat.fittnamefixp == name)[0]
            factfixpplot = gdat.fittfactfixpplot[indxfixptemp]
            setattr(gdat, 'fact' + name + 'plot', factfixpplot)
        else:
            try:    
                getattr(gdat, 'fact' + name + 'plot')
            except:
                setattr(gdat, 'fact' + name + 'plot', 1.)

    # copy fixp variables to individual variables
    for k, namefixp in enumerate(gdat.fittnamefixp):
        setattr(gdat, 'fittminm' + namefixp, gdat.fittminmfixp[k])
        setattr(gdat, 'fittmaxm' + namefixp, gdat.fittmaxmfixp[k])
        setattr(gdat, 'fittscal' + namefixp, gdat.fittscalfixp[k])
    
    # log-prior register
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    numb = 0
    for l in gdat.fittindxpopl:
        numb += len(gdat.fittliststrgfeat[l])
    
    # process index
    gdat.indxproc = np.arange(gdat.numbproc)

    if gdat.fittlensmodltype != 'none' or gdat.datatype == 'mock' and gdat.truelensmodltype != 'none':
        retr_axis(gdat, 'mcut', gdat.minmmcut, gdat.maxmmcut, gdat.numbbinsplot)
        retr_axis(gdat, 'bein', gdat.minmbein, gdat.maxmbein, gdat.numbbinsplot)

    # angular deviation
    gdat.numbanglhalf = 10
    gdat.indxanglhalf = np.arange(gdat.numbanglhalf)
    retr_axis(gdat, 'anglhalf', 0., gdat.maxmgangdata, gdat.numbanglhalf)
    gdat.numbanglfull = 1000
    gdat.indxanglfull = np.arange(gdat.numbanglfull)
    retr_axis(gdat, 'anglfull', 0., 3. * gdat.maxmgangdata, gdat.numbanglfull)
    
    # temp
    #gdat.binsanglcosi = np.sort(np.cos(gdat.binsangl))
    
    # temp
    #gdat.meshbackener = np.meshgrid(gdat.indxback, gdat.indxener, indexing='ij')
    
    # plotting
    ## the normalized offset for text annotation of point sources in the frames
    gdat.offstextimag = gdat.maxmgangdata * 0.05
    
    ## figure size
    gdat.plotsize = 12
    ## size of the images
    gdat.sizeimag = 1.3 * gdat.plotsize
    ## text
    gdat.fittlegd = 'Model'
    if gdat.datatype == 'mock':
        gdat.refrlegd = 'True'
    else:
        gdat.refrlegd = 'Ref'
    for strgmodl in ['refr', 'fitt']:
        
        legdpopl = getattr(gdat, strgmodl + 'legdpopl')
        
        legd = getattr(gdat, strgmodl + 'legd')
        if strgmodl == 'refr':
            indxpopl = gdat.indxrefr
        else:
            indxpopl = gdat.fittindxpopl

        if legdpopl is None:
            legdpopl = [[] for l in indxpopl]
            for l in indxpopl:
                legdpopl[l] = 'Pop %d' % l
        
        legdelem = [[] for l in indxpopl]
        legdmiss = [[] for l in indxpopl]
        legdhits = [[] for l in indxpopl]
        for l in indxpopl:
            legdelem[l] = legd + ' ' + legdpopl[l]
            legdmiss[l] = legdelem[l] + ' miss'
            legdhits[l] = legdelem[l] + ' hit'
        legdhost = legd + ' host'
        legdsour = legd + ' sour'
        setp_varbvalu(gdat, 'legdelem', legdelem, strgmodl=strgmodl)
        setattr(gdat, strgmodl + 'legdmiss', legdmiss)
        setattr(gdat, strgmodl + 'legdhits', legdhits)
        setattr(gdat, strgmodl + 'legdhost', legdhost)
        setattr(gdat, strgmodl + 'legdsour', legdsour)

    ## PSF class indices for which images will be plotted
    if gdat.numbevtt == 1:
        gdat.indxevttplot = gdat.indxevtt
    else:
        gdat.indxevttplot = np.concatenate((np.array([-1]), gdat.indxevtt))
    
    gdat.numbenerevtt = gdat.numbener * gdat.numbevtt
    
    gdat.listnamechro = ['totl', 'prop', 'diag', 'save', 'plot', 'proc', 'pars', 'modl', 'llik', 'sbrtmodl']
    gdat.listlegdchro = ['Total', 'Proposal', 'Diagnostics', 'Save', 'Plot', 'Process', 'Parse', 'Model', 'Likelihood', 'Total emission']
    if gdat.fittnumbtrap > 0:
        gdat.listnamechro += ['evalelem']
        gdat.listlegdchro += ['El. Eval. Arrays']
    if gdat.fittlensmodltype != 'none':
        gdat.listnamechro += ['deflzero', 'deflhost', 'deflextr', 'sbrtlens', 'sbrthost']
        gdat.listlegdchro += ['Array initialization', 'Host Deflection', 'External deflection', 'Lensed emission', 'Host emission']
    if gdat.fittboolelemsbrtdfncanyy:
        gdat.listnamechro += ['elemsbrtdfnc']
        gdat.listlegdchro += ['Dfnc S Brght']
    if gdat.fittboolelemdeflsubhanyy:
        gdat.listnamechro += ['elemdeflsubh']
        gdat.listlegdchro += ['Subh Defl']
    if gdat.fittboolelemsbrtextsbgrdanyy:
        gdat.listnamechro += ['elemsbrtextsbgrd']
        gdat.listlegdchro += ['Bkg Exts S Brght']
    booltemp = False
    for strgmodl in gdat.liststrgmodl:
        booltemp = booltemp or getattr(gdat, strgmodl + 'psfnevaltype')
    if booltemp or gdat.fittpsfnevaltype == 'full' or gdat.truepsfnevaltype == 'full':
        gdat.listnamechro += ['psfnconv']
        gdat.listlegdchro += ['Img for PSF Conv.']
    
    gdat.listnamechro += ['expo', 'lpri', 'tert']
    gdat.listlegdchro += ['Exposure', 'Prior', 'Tertiary']
    gdat.numbchro = len(gdat.listnamechro)
    
    # temp
    gdat.boolintpanglcosi = False

    if gdat.boolthindata:
        gdat.factdatathin = 10
        if gdat.pixltype != 'cart' or gdat.numbsidecart % gdat.factdatathin != 0:
            raise Exception('Cannot thin the data.')
        #gdat.indxpixlkeep = gdat.indxpixlfull[::gdat.factdatathin]
        #gdat.numbpixlkeep = gdat.indxpixlkeep.size
        gdat.indxpixlkill = np.setdiff1d(gdat.indxpixlfull, gdat.indxpixlkeep)
        gdat.numbsidecart = gdat.numbsidecart / 10
        gdat.lgalgrid = gdat.lgalgrid[gdat.indxpixlkeep]
        gdat.bgalgrid = gdat.bgalgrid[gdat.indxpixlkeep]
        gdat.indxpixlfull = gdat.indxpixlfull[gdat.indxpixlkeep]
        
    # the function to measure time
    # temp
    gdat.strgfunctime = 'clck'
    if gdat.strgfunctime == 'clck':
        gdat.functime = time.clock
    if gdat.strgfunctime == 'time':
        gdat.functime = time.time

    ## longitude
    gdat.numblgalpntsprob = gdat.numbsidepntsprob
    gdat.numbbgalpntsprob = gdat.numbsidepntsprob
    gdat.binslgalpntsprob = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.binsbgalpntsprob = np.linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.indxlgalpntsprob = np.arange(gdat.numblgalpntsprob)
    gdat.indxbgalpntsprob = np.arange(gdat.numbbgalpntsprob)

    retr_axis(gdat, 'defl', -gdat.maxmgangdata, gdat.maxmgangdata, 40)
    retr_axis(gdat, 'deflsubh', -gdat.maxmgangdata * 1e-2, gdat.maxmgangdata * 1e-2, 40)

    # lensing problem setup
    ## number of deflection components to plot

    gdat.binslgalcartmesh, gdat.binsbgalcartmesh = np.meshgrid(gdat.binslgalcart, gdat.binsbgalcart, indexing='ij')
    gdat.meanlgalcartmesh, gdat.meanbgalcartmesh = np.meshgrid(gdat.meanlgalcart, gdat.meanbgalcart, indexing='ij')
    if gdat.pixltype == 'cart':
        gdat.sizepixl = np.sqrt(gdat.apix)
        gdat.indxsidecart = np.arange(gdat.numbsidecart)
        gdat.indxpixlrofi = np.arange(gdat.numbpixlcart)
        gdat.indxsidemesh = np.meshgrid(gdat.indxsidecart, gdat.indxsidecart, indexing='ij')
        gdat.lgalgrid = gdat.meanlgalcart[gdat.indxsidemesh[0].flatten()]
        gdat.bgalgrid = gdat.meanbgalcart[gdat.indxsidemesh[1].flatten()]
        gdat.shapcart = (gdat.numbsidecart, gdat.numbsidecart)
        gdat.lgalgridfull = np.copy(gdat.lgalgrid)
        gdat.bgalgridfull = np.copy(gdat.bgalgrid)
        gdat.lgalgridcart = gdat.lgalgrid.reshape(gdat.shapcart)
        gdat.bgalgridcart = gdat.bgalgrid.reshape(gdat.shapcart)
        gdat.indxpent = np.meshgrid(gdat.indxener, gdat.indxsidecart, gdat.indxsidecart, gdat.indxevtt, indexing='ij')
    if gdat.pixltype == 'heal':
        lgalheal, bgalheal, gdat.numbpixlfull, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)
        lgalheal = np.deg2rad(lgalheal)
        bgalheal = np.deg2rad(bgalheal)
   
        gdat.indxpixlrofi = np.where((np.fabs(lgalheal) < gdat.maxmgangdata) & (np.fabs(bgalheal) < gdat.maxmgangdata))[0]
        
        gdat.indxpixlrofimarg = np.where((np.fabs(lgalheal) < 1.2 * gdat.maxmgangdata) & (np.fabs(bgalheal) < 1.2 * gdat.maxmgangdata))[0]

        gdat.lgalgrid = lgalheal
        gdat.bgalgrid = bgalheal
    
    gdat.indxpixlfull = np.arange(gdat.numbpixlfull)
    if gdat.pixltype == 'cart':
        gdat.indxpixlcart = np.arange(gdat.numbpixlcart)
    
    if gdat.evttbins:
        # PSF class string
        gdat.strgevtt = []
        for m in gdat.indxevtt:
            gdat.strgevtt.append('PSF%d' % gdat.indxevttincl[m])
    
    # angular diameter distance
    #if gdat.fittlensmodltype != 'none' or gdat.datatype == 'mock' and gdat.truelensmodltype != 'none':
    #    gdat.adis = gdat.adishost
    
    # power spectra
    if gdat.pixltype == 'cart':
        gdat.numbwvecodim = gdat.numbsidecart
        gdat.minmanglodim = 0.
        gdat.maxmanglodim = 2. * gdat.maxmgangdata
        gdat.minmmpolodim = 0.
        gdat.maxmmpolodim = 1. / 2. / gdat.sizepixl
        retr_axis(gdat, 'anglodim', gdat.minmanglodim, gdat.maxmanglodim, gdat.numbsidecart, invr=True)
        retr_axis(gdat, 'mpolodim', gdat.minmmpolodim, gdat.maxmmpolodim, gdat.numbsidecart / 2)
        if gdat.fittlensmodltype != 'none' or gdat.datatype == 'mock' and gdat.truelensmodltype != 'none':
            # temp -- this should minima, maxima of adishost and the true metamodel into account
            gdat.minmwvecodim = gdat.minmmpolodim / np.amax(gdat.fittadishost)
            gdat.maxmwvecodim = gdat.maxmmpolodim / np.amin(gdat.fittadishost)
            gdat.minmwlenodim = gdat.minmanglodim * np.amin(gdat.fittadishost)
            gdat.maxmwlenodim = gdat.maxmanglodim * np.amax(gdat.fittadishost)
            retr_axis(gdat, 'wvecodim', gdat.minmwvecodim, gdat.maxmwvecodim, gdat.numbsidecart / 2)
            retr_axis(gdat, 'wlenodim', gdat.minmwlenodim, gdat.maxmwlenodim, gdat.numbsidecart, invr=True)
            gdat.meanwveclgal, gdat.meanwvecbgal = np.meshgrid(gdat.meanwvecodim, gdat.meanwvecodim, indexing='ij')
            gdat.meanwvec = np.sqrt(gdat.meanwveclgal**2 + gdat.meanwvecbgal**2)
        gdat.meanmpollgal, gdat.meanmpolbgal = np.meshgrid(gdat.meanmpolodim, gdat.meanmpolodim, indexing='ij')
        gdat.meanmpol = np.sqrt(gdat.meanmpollgal**2 + gdat.meanmpolbgal**2)

    # element parameter vector indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompflux = 2
    gdat.indxcompsind = 3
    gdat.indxcompcurv = 4
    gdat.indxcompexpc = 4

    # check the exposure map data structure
    booltemp = False
    if gdat.expo.ndim != 3:
        booltemp = True
    if gdat.pixltype == 'cart' and gdat.expo.shape[1] != gdat.numbpixlcart:
        booltemp = True
    if booltemp:
        print 'gdat.expo'
        summgene(gdat.expo)
        print 'gdat.numbsidecart'
        print gdat.numbsidecart
        raise Exception('Exposure does not have the right data structure. It should be a list of 3D np.arrays.')

    if gdat.boolsqzeexpo:
        gdat.expo *= 1e-10
    if gdat.boolexplexpo:
        gdat.expo *= 1e10
    
    gdat.factcmplplot = 1.
    gdat.factfdisplot = 1.

    if gdat.boolthindata:
        #gdat.expo[:, gdat.indxpixlkill, :] = 0.
        expotemp = np.copy(gdat.expo[:, gdat.indxpixlfull[::gdat.factdatathin], :])
        sbrttemp = np.copy(gdat.sbrtdata[:, gdat.indxpixlfull[::gdat.factdatathin], :])
        gdat.expo = expotemp 
        gdat.sbrtdata = sbrttemp
        
    # only include desired energy and PSF class bins
    gdat.indxcubeincl = np.meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    
    ## exposure
    if gdat.correxpo:
        # temp -- for some reason lists of np.arrays require manual processing
        gdat.expo = gdat.expo[gdat.indxcubeincl]
        if gdat.datatype == 'inpt':
            gdat.sbrtdata = gdat.sbrtdata[gdat.indxcubeincl]
    
    ## backgrounds
    for strgmodl in gdat.liststrgmodl:
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        indxback = getattr(gdat, strgmodl + 'indxback')
        sbrtbacknormincl = [[] for c in indxback]
        for c in indxback:
            sbrtbacknormincl[c] = sbrtbacknorm[c][gdat.indxcubeincl]
        setattr(gdat, strgmodl + 'sbrtbacknorm', sbrtbacknormincl)
    
    # obtain cartesian versions of the maps
    #if gdat.pixltype == 'cart':
    #    gdat.expocart = gdat.expo.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
    #    for strgmodl in gdat.liststrgmodl:
    #        sbrtbacknormcart = []
    #        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
    #        for c in getattr(gdat, strgmodl + 'indxback'):
    #            sbrtbacknormcart.append(sbrtbacknorm[c].reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt)))
    #        setattr(gdat, strgmodl + 'sbrtbacknormcart', sbrtbacknormcart)
    
    # mask the exposure map
    if gdat.listmask is not None:
        for mask in gdat.listmask:
            if mask[0] == 'sqre':
                indxpixlmask = np.where((gdat.lgalgrid > mask[1]) & (gdat.lgalgrid < mask[2]) & (gdat.bgalgrid > mask[3]) & (gdat.bgalgrid < mask[4]))[0]
            if mask[0] == 'circ':
                indxpixlmask = np.where(np.sqrt((gdat.lgalgrid - mask[1])**2 + (gdat.bgalgrid - mask[2])**2) < mask[3])[0]
            if mask[0] == 'hstr':
                indxpixlmask = np.where((gdat.bgalgrid > mask[1]) & (gdat.bgalgrid < mask[2]))[0]
            if gdat.masktype == 'zero':
                gdat.expo[:, indxpixlmask, :] = 0.
            if gdat.masktype == 'ignr':
                gdat.expo[:, indxpixlmask, :] = 1e-49

    # plotting
    ## ROI
    if gdat.numbpixlfull != 1:
        gdat.exttrofi = np.array([gdat.minmlgaldata, gdat.maxmlgaldata, gdat.minmbgaldata, gdat.maxmbgaldata])
        gdat.exttrofi *= gdat.anglfact 
        gdat.frambndrdata = gdat.maxmgangdata * gdat.anglfact

    ## marker size
    gdat.minmmrkrsize = 100
    gdat.maxmmrkrsize = 500
    ## marker line width
    gdat.mrkrlinewdth = 3
    ## marker opacity
    gdat.alphhist = 0.5
    gdat.alphline = 0.5
    gdat.alphbndr = 0.5
    gdat.alphelem = 1.
    gdat.alphmaps = 1.
    
    # number of colorbar ticks in the maps
    gdat.numbtickcbar = 11
    
    ## color bars
    gdat.minmlpdfspatpriointp = np.log(1. / 2. / gdat.maxmgangdata) - 10.
    gdat.maxmlpdfspatpriointp = np.log(1. / 2. / gdat.maxmgangdata) + 10.
    gdat.scallpdfspatpriointp = 'self'
    gdat.cmaplpdfspatpriointp = 'PuBu'
    
    gdat.minmllikmaps = -10.
    gdat.maxmllikmaps = 0.
    gdat.scalllikmaps = 'asnh'
    gdat.cmapllikmaps = 'YlGn'
    
    gdat.minmperc = 0.
    gdat.maxmperc = 1e2
    gdat.scalperc = 'asnh'
    gdat.cmapperc = 'afmhot'
    
    gdat.minmpercresi = -1e2
    gdat.maxmpercresi = 1e2
    gdat.scalpercresi = 'asnh'
    gdat.cmappercresi = 'coolwarm'
    
    gdat.scalexpo = 'logt'
    gdat.cmapexpo = 'OrRd'
    
    gdat.scalcntpdata = 'logt'
    gdat.cmapcntpdata = 'Greys'
    
    gdat.scalcntpmodl = 'logt'
    gdat.cmapcntpmodl = 'Greys'
    
    gdat.scalcntpresi = 'asnh'
    gdat.cmapcntpresi = make_cmapdivg('Red', 'Orange')

    gdat.minmconv = 1e-2
    gdat.maxmconv = 10.
    gdat.scalconv = 'logt'
    gdat.cmapconv = 'Purples'
    
    gdat.minmconvelem = 1e-4
    gdat.maxmconvelem = 1e-1
    gdat.scalconvelem = 'logt'
    gdat.cmapconvelem = 'Purples'
    
    gdat.minms2nr = 0.
    gdat.maxms2nr = 10.
    gdat.scals2nr = 'asnh'
    gdat.cmaps2nr = 'magma'
    
    gdat.minmmagn = -1e2
    gdat.maxmmagn = 1e2
    gdat.scalmagn = 'asnh'
    gdat.cmapmagn = 'BrBG'
    
    gdat.minmdeflresiperc = -100.
    gdat.maxmdeflresiperc = 100.
    gdat.scaldeflresiperc = 'self'
    gdat.cmapdeflresiperc = 'Oranges'
    
    gdat.minmconvelemresi = -0.1
    gdat.maxmconvelemresi = 0.1
    gdat.scalconvelemresi = 'self'
    gdat.cmapconvelemresi = 'PiYG'
    
    gdat.minmconvelemresiperc = -100.
    gdat.maxmconvelemresiperc = 100.
    gdat.scalconvelemresiperc = 'self'
    gdat.cmapconvelemresiperc = 'PiYG'
    
    gdat.minmmagnresi = -10.
    gdat.maxmmagnresi = 10.
    gdat.scalmagnresi = 'self'
    gdat.cmapmagnresi = 'PRGn'
    
    gdat.minmmagnresiperc = -100.
    gdat.maxmmagnresiperc = 100.
    gdat.scalmagnresiperc = 'self'
    gdat.cmapmagnresiperc = 'PRGn'
    
    gdat.lgalgrid = gdat.lgalgrid[gdat.indxpixlrofi]
    gdat.bgalgrid = gdat.bgalgrid[gdat.indxpixlrofi]
   
    if np.amax(gdat.expo) <= 0.:
        raise Exception('Bad exposure.')

    # temp
    #gdat.expo[np.where(gdat.expo < 1e-50)] = 1e-50
    
    # exclude voxels with vanishing exposure
    if gdat.correxpo:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.indxpixlrofi = np.intersect1d(gdat.indxpixlrofi, np.where(gdat.expo[i, :, m] > 0.)[0])
    
    gdat.indxcuberofi = np.meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = np.arange(gdat.numbpixl)
    gdat.numbdata = gdat.numbener * gdat.numbevtt * gdat.numbpixl

    #gdat.lgalgridrofi = gdat.lgalgrid[gdat.indxpixlrofi]
    #gdat.bgalgridrofi = gdat.bgalgrid[gdat.indxpixlrofi]

    gdat.minmexpo = np.amin(gdat.expo[np.where(gdat.expo > 1e-100)])
    gdat.maxmexpo = np.amax(gdat.expo)
    gdat.minmexpo = np.amin(gdat.minmexpo)
    gdat.maxmexpo = np.amax(gdat.maxmexpo)

    if gdat.minmexpo > 0:
        gdat.indxpixlroficnvt = np.arange(gdat.numbpixlfull)
    else:
        if gdat.verbtype > 0:
            print 'Calculating lookup table for zero-exposure pixels...'
        cntr = 0
        gdat.indxpixlroficnvt = full(gdat.numbpixlfull, -1)
        for j in gdat.indxpixlfull:
            if j in gdat.indxpixlrofi:
                gdat.indxpixlroficnvt[j] = cntr
                cntr += 1
    
    if gdat.datatype == 'inpt':
        gdat.sbrtdata = gdat.sbrtdata[gdat.indxcuberofi]

    ## exposure
    if gdat.correxpo:
        gdat.expofull = np.copy(gdat.expo)
        gdat.expo = gdat.expo[gdat.indxcuberofi]
    
    ## backgrounds
    for strgmodl in gdat.liststrgmodl:
        indxback = getattr(gdat, strgmodl + 'indxback')
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        if gdat.pixltype == 'heal':
            sbrtbackhealfull = [[] for c in indxback]
            for c in indxback:
                sbrtbackhealfull[c] = np.copy(sbrtbacknorm[c])
            setattr(gdat, strgmodl + 'sbrtbackhealfull', sbrtbackhealfull)
        sbrtbacknormincl = [[] for c in indxback]
        for c in indxback:
            sbrtbacknormincl[c] = sbrtbacknorm[c][gdat.indxcuberofi]
        setattr(gdat, strgmodl + 'sbrtbacknorm', sbrtbacknormincl)
                
    gdat.expototl = []
    gdat.expototlmean = []
    gdat.expototl = np.sum(gdat.expo, axis=2)
    gdat.expototlmean = np.mean(gdat.expototl, axis=1)

    if 'locl' in gdat.commelemspatevaltype:
        if gdat.exprtype == 'sdyn':
            gdat.maxmangl = 1.
        if gdat.exprtype == 'ferm':
            gdat.maxmangl = 20. / gdat.anglfact
        if gdat.exprtype == 'tess':
            gdat.maxmangl = 25. / gdat.anglfact
        if gdat.exprtype == 'chan':
            gdat.maxmangl = 15. / gdat.anglfact
        if gdat.exprtype == 'hubb':
            gdat.maxmangl = 1. / gdat.anglfact
    else:
        gdat.maxmangl = gdat.maxmgangdata * np.sqrt(2.) * 2. * 1.1
        
    if gdat.numbpixlfull != 1:
        retr_axis(gdat, 'angl', 0., gdat.maxmangl, gdat.numbangl)
        
    gdat.listnamespatmean = ['full']
    if gdat.exprtype == 'ferm':
        gdat.listnamespatmean += ['innr']
    gdat.numbspatmean = len(gdat.listnamespatmean)
    gdat.indxspatmean = np.arange(gdat.numbspatmean)
    gdat.listindxcubespatmean = [[] for b in gdat.indxspatmean]
    gdat.indxcube = np.meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        if namespatmean == 'full':
            gdat.listindxcubespatmean[b] = gdat.indxcube
        if namespatmean == 'innr':
            gdat.indxpixlinnr = np.where(np.sqrt(gdat.lgalgrid**2 + gdat.bgalgrid**2) < 5. / gdat.anglfact)[0]
            gdat.listindxcubespatmean[b] = np.meshgrid(gdat.indxener, gdat.indxpixlinnr, gdat.indxevtt, indexing='ij')
    
    if gdat.numbpixl > 1:
        # store pixels as unit vectors
        gdat.xdatgrid, gdat.ydatgrid, gdat.zaxigrid = retr_unit(gdat.lgalgrid, gdat.bgalgrid)
   
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
                gdat.pixlcnvt = np.zeros(gdat.numbpixlfull, dtype=int) - 1
                numbpixlmarg = gdat.indxpixlrofimarg.size
                for k in range(numbpixlmarg):
                    dist = retr_angldistunit(gdat, lgalheal[gdat.indxpixlrofimarg[k]], bgalheal[gdat.indxpixlrofimarg[k]], gdat.indxpixl)
                    gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)
                fobj = open(path, 'wb')
                cPickle.dump(gdat.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
                fobj.close()
        
        # dummy pixel indices for full (nonlocal) element kernel evaluation 
        gdat.listindxpixl = []
        if gdat.datatype == 'mock':
            numb = max(np.sum(gdat.fittmaxmnumbelem), np.sum(gdat.truemaxmnumbelem)) + 2
        else:
            numb = np.sum(gdat.fittmaxmnumbelem) + 2
        for k in range(int(numb)):
            gdat.listindxpixl.append([])
            for kk in range(k):
                gdat.listindxpixl[k].append(gdat.indxpixl)
        
        # spatial averaging setup
        # temp
        
        # temp -- check if 1000 is too much
        gdat.numbanglelem = 1000
    
    # turn off relevant proposal types
    gdat.numbprop = 5
    gdat.indxprop = np.arange(gdat.numbprop)
    indxelemfull = [range(gdat.fittmaxmnumbelem[l]) for l in gdat.fittindxpopl]
    indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, 'fitt')
    
    gdat.numbstdp = gdat.fittnumbfixp - gdat.fittnumbpopl
    cntr = 0
    for l in gdat.fittindxpopl:
        for strgcomp in gdat.fittliststrgcomp[l]:
            setattr(gdat, 'indxstdppop%d' % l + strgcomp, gdat.numbstdp)
            cntr += 1
    gdat.numbstdp += cntr
    
    gdat.strgstdp = np.copy(np.array(gdat.fittlablfixp[gdat.fittnumbpopl:]))
    gdat.namestdp = np.copy(np.array(gdat.fittnamefixp[gdat.fittnumbpopl:]))
    for l in gdat.fittindxpopl:
        for strgcomp in gdat.fittliststrgcomp[l]:
            gdat.strgstdp = np.append(gdat.strgstdp, getattr(gdat, 'labl' + strgcomp))
            gdat.namestdp = np.append(gdat.namestdp, strgcomp + 'pop%d' % l)
    gdat.strgstdp = list(gdat.strgstdp)
    gdat.indxstdp = np.arange(gdat.numbstdp)
    gdat.indxstdpprop = gdat.indxstdp
    
    # proposal scale indices for each parameter
    gdat.indxstdppara = np.zeros(gdat.fittnumbpara, dtype=int) - 1
    cntr = 0
    gdat.indxstdppara[gdat.fittnumbpopl:gdat.fittnumbfixp] = gdat.fittindxfixp[gdat.fittnumbpopl:] - gdat.fittnumbpopl
    if gdat.fittnumbtrap > 0:
        for l in gdat.fittindxpopl:
            for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                for indx in indxsampcomp[strgcomp][l]:
                    gdat.indxstdppara[indx] = cntr + gdat.fittnumbfixp - gdat.fittnumbpopl
                cntr += 1
    
    # for the fitting model, define proposal type indices
    dicttemptemp = deepcopy(gdat.__dict__) 
    for name, valu in dicttemptemp.iteritems():
        if name.startswith('fittindxfixp') and name != 'fittindxfixp' and not name.startswith('fittindxfixpnumbelem'):
            indxstdp = gdat.indxstdppara[valu]
            setattr(gdat, 'indxstdp' + name[12:], indxstdp)
    
    # proposal scale
    if gdat.fittlensmodltype != 'none' or gdat.datatype == 'mock' and gdat.truelensmodltype != 'none':
        
        gdat.stdp = 1e-4 + np.zeros(gdat.numbstdp)
        
        if gdat.fittnumbpopl > 0:
            if gdat.fittmaxmnumbelempopl[0] > 0:
                gdat.stdp[gdat.indxstdpmeanelempop0] = 1e-1

        gdat.stdp[gdat.indxstdpsigcen00evt0] = 3e-2
        gdat.stdp[gdat.indxstdpbacpback0000en00] = 1e-3
        gdat.stdp[gdat.indxstdpbacpback0000en00] = 1e-1
        
        if gdat.fittlensmodltype != 'none':
            gdat.stdp[gdat.indxstdplgalsour] = 1e-3
            gdat.stdp[gdat.indxstdpbgalsour] = 1e-3
            gdat.stdp[gdat.indxstdpfluxsour] = 1e-2
            if gdat.numbener > 1:
                gdat.stdp[gdat.indxstdpsindsour] = 1e-3
            gdat.stdp[gdat.indxstdpsizesour] = 1e-1
            gdat.stdp[gdat.indxstdpellpsour] = 1e-1
            gdat.stdp[gdat.indxstdpanglsour] = 1e-1
        if gdat.fitthostemistype != 'none':
            gdat.stdp[gdat.indxstdplgalhostisf0] = 3e-4
            gdat.stdp[gdat.indxstdpbgalhostisf0] = 3e-4
            gdat.stdp[gdat.indxstdpfluxhostisf0] = 1e-3
            if gdat.numbener > 1:
                gdat.stdp[gdat.indxstdpsindhostisf0] = 1e-3
            gdat.stdp[gdat.indxstdpsizehostisf0] = 3e-3
        if gdat.fittlensmodltype != 'none':
            gdat.stdp[gdat.indxstdpbeinhostisf0] = 1e-3
        if gdat.fitthostemistype != 'none':
            gdat.stdp[gdat.indxstdpellphostisf0] = 1e-2
            gdat.stdp[gdat.indxstdpanglhostisf0] = 1e-2
            gdat.stdp[gdat.indxstdpserihostisf0] = 1e-2
        if gdat.fittlensmodltype != 'none':
            gdat.stdp[gdat.indxstdpsherextr] = 1e-1
            gdat.stdp[gdat.indxstdpsangextr] = 3e-2
        
    else:
        
        if gdat.exprtype == 'ferm':
            gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
            
            if gdat.fittnumbtrap > 0:
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpmeanelem]] = 4e-2
                
                for l in gdat.fittindxpopl:
                    if gdat.fittfluxdisttype[l] == 'powrslop':
                        gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpfluxdistsloppop0]] = 1e-1
                    else:
                        gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpfluxdistbrekpop0]] = 1e-1
                        gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpfluxdistsloplowrpop0]] = 1e-1
                        gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpfluxdistslopupprpop0]] = 1e-1
            
            gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en00')]] = 5e-3
            gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en01')]] = 1e-2
            gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en02')]] = 3e-2
            
            if 'fdfm' in gdat.fittlistnameback: 
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0001en00')]] = 8e-4
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0001en01')]] = 1e-3
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0001en02')]] = 2e-3
            
            if 'dark' in gdat.fittlistnameback: 
                indxbackdark = gdat.fittlistnameback.index('dark')
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback%04den00' % indxbackdark)]] = 2e-2
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback%04den01' % indxbackdark)]] = 2e-2
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback%04den02' % indxbackdark)]] = 3e-2
            
            if gdat.fittnumbtrap > 0:
                if 'excefixd' in gdat.strgcnfg:
                    gdat.stdp[gdat.indxstdppop0lgal] = 1e-100
                    gdat.stdp[gdat.indxstdppop0bgal] = 1e-100
                gdat.stdp[gdat.indxstdppop0flux] = 8e-2

                if gdat.fittspectype[0] == 'colr':
                    gdat.stdp[gdat.indxstdppop0sindcolr0001] = 8e-2
                    gdat.stdp[gdat.indxstdppop0sindcolr0002] = 2e-1
        
        if gdat.exprtype == 'chan':
            gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
            
            if gdat.fittnumbtrap > 0:
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpmeanelem]] = 2e-1
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpfluxdistsloppop0]] = 2e-1
            
            if gdat.fittnumbtrap > 0 and gdat.numbpixlfull > 1:
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixppsfp]] = 4e-1
            
            if gdat.indxenerincl.size == 5 and (gdat.indxenerincl == np.arange(5)).all():
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en00')]] = 2e-2
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en01')]] = 3e-2
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en02')]] = 2e-2
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en03')]] = 2e-2
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en04')]] = 1e-2
            elif gdat.indxenerincl.size == 2 and (gdat.indxenerincl == np.array([2])).all():
                gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en00')]] = 2e-2
            
            if gdat.fittnumbtrap > 0:
                if gdat.numbpixlfull > 1:
                    gdat.stdp[gdat.indxstdppop0lgal] = 2e-2
                    gdat.stdp[gdat.indxstdppop0bgal] = 2e-2
                    if gdat.numbener > 1:
                        gdat.stdp[gdat.indxstdppop0sind] = 2e-1
                gdat.stdp[gdat.indxstdppop0flux] = 2e-1
        
        if gdat.exprtype == 'sdyn':
            gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
        
            if gdat.fittnumbtrap > 0:
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpmeanelem]] = 2e-1
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpnobjdistsloppop0]] = 3e-1
                try:
                    gdat.stdp[gdat.indxstdppara[gdat.fittindxfixpgwdtdistsloppop0]] = 3e-1
                except:
                    pass

            if gdat.fittpsfnevaltype != 'none':
                gdat.stdp[gdat.indxstdppara[gdat.fittindxfixppsfp]] = 4e-1
            
            gdat.stdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback0000en00')]] = 2e-2
        
            if gdat.fittnumbtrap > 0:
                gdat.stdp[gdat.indxstdppop0lgal] = 4e-2
                gdat.stdp[gdat.indxstdppop0bgal] = 4e-2
                gdat.stdp[gdat.indxstdppop0nobj] = 3e-1
                try:
                    gdat.stdp[gdat.indxstdppop0gwdt] = 5e-1
                except:
                    pass

        if gdat.exprtype == 'fire':
            gdat.stdp = 1e-2 + np.zeros(gdat.numbstdp)
    
    if gdat.boolsqzeprop:
        gdat.stdp[:]= 1e-100
    
    if gdat.boolexplprop:
        gdat.stdp[:] = 1e100

    if (gdat.stdp > 1e100).any():
        raise Exception('')
        
    if (gdat.stdp == 0).any():
        raise Exception('')
        
    if gdat.stdp.size != gdat.numbstdp or gdat.indxstdp.size != gdat.stdp.size:
        print 'gdat.numbstdp'
        print gdat.numbstdp
        print 'gdat.stdp'
        print gdat.stdp
        summgene(gdatmodi.stdp)
        raise Exception('')

    if gdat.verbtype > 1 and boolinitsetp:
        # temp
        for strgmodl in gdat.liststrgmodl:
            print 'strgmodl'
            print strgmodl
            print 'Fixed dimensional parameters:'
            print '%20s%25s%5s%20s%20s' % ('name', 'labl', 'scal', 'minm', 'maxm')
            for k in getattr(gdat, strgmodl + 'indxfixp'):
                namefixp = getattr(gdat, strgmodl + 'namefixp')[k]
                lablfixp = getattr(gdat, strgmodl + 'lablfixp')[k]
                scalfixp = getattr(gdat, strgmodl + 'scalfixp')[k]
                minmfixp = getattr(gdat, strgmodl + 'minmfixp')[k]
                maxmfixp = getattr(gdat, strgmodl + 'maxmfixp')[k]
                print '%20s%25s%5s%20.6g%20.6g' % (namefixp, lablfixp, scalfixp, minmfixp, maxmfixp)
            
            indxpopl = getattr(gdat, strgmodl + 'indxpopl')
            liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
            listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
            liststrgfeatmodu = getattr(gdat, strgmodl + 'liststrgfeatmodu')
            liststrgpdfnmodu = getattr(gdat, strgmodl + 'liststrgpdfnmodu')
            liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
            liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
            
            print 'Element parameters'
            print '%20s%20s' % ('strgcomp', 'scalcomp')
            for l in indxpopl:
                for strgcomp, scalcomp in zip(liststrgcomp[l], listscalcomp[l]):
                    print '%20s%20s' % (strgcomp, scalcomp)
            
            print '%20s%20s' % ('strgmodu', 'pdfnmodu')
            for l in indxpopl:
                for strgmodu, pdfnmodu in zip(liststrgfeatmodu[l], liststrgpdfnmodu[l]):
                    print '%20s%20s' % (strgmodu, pdfnmodu)
            
            print '%20s%20s' % ('strgfeat', 'pdfnprio')
            for l in indxpopl:
                for strgfeat, pdfnprio in zip(liststrgfeatprio[l], liststrgpdfnprio[l]):
                    print '%20s%20s' % (strgfeat, pdfnprio)
            
    # proposals
    # parameters not subject to proposals
    
    if gdat.probtran is None:
        if gdat.fittnumbtrap > 0:
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
    if gdat.probspmr is None:
        if gdat.fittnumbtrap > 0:
            gdat.probspmr = gdat.probtran / 2.
        else:
            gdat.probspmr = 0.
    
    gdat.probbrde = 1. - gdat.probspmr

    if gdat.probbrde < 0:
        raise Exception('')
    gdat.legdproptype = ['Within']
    gdat.nameproptype = ['with']
    print 'gdat.fittnumbtrap'
    print gdat.fittnumbtrap
    if gdat.fittnumbtrap > 0:
        gdat.legdproptype += ['Birth', 'Death', 'Split', 'Merge']
        gdat.nameproptype += ['brth', 'deth', 'splt', 'merg']
    gdat.numbproptype = len(gdat.legdproptype)
    gdat.nameproptype = np.array(gdat.nameproptype)
    cntr = tdpy.util.cntr()
    if gdat.fittnumbtrap > 0.:
        # birth
        gdat.indxproptypebrth = cntr.incr()
        # death
        gdat.indxproptypedeth = cntr.incr()
        if gdat.probspmr > 0.:
            # split
            gdat.indxproptypesplt = cntr.incr()
            # merge
            gdat.indxproptypemerg = cntr.incr()
   
    gdat.indxproptype = np.arange(gdat.numbproptype)
    gdat.indxfixpprop = np.arange(gdat.fittnumbfixp)
    gdat.numbstdpfixp = gdat.fittnumbfixp - gdat.fittnumbpopl
    #### filter for model elements
    gdat.listnamefilt = ['']
    if gdat.priofactdoff != 1.:
        gdat.listnamefilt += ['pars']
    #### model elements inside the image
    if gdat.commboolelempsfnanyy:
        gdat.listnamefilt += ['bndr']
    #### model subhalos inside high normalized relevance region
    if 'lens' in gdat.commelemtype:
        gdat.listnamefilt += ['nrel']
    
    if gdat.datatype == 'inpt':
        proc_cntpdata(gdat)
    
    # set common plotting minima and maxima for element features
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
        for l in indxpopl:
            for namefeat in liststrgfeat[l]:
                try:
                    getattr(gdat, 'minm' + namefeat)
                    getattr(gdat, 'maxm' + namefeat)
                except:
                    try:
                        setattr(gdat, 'minm' + namefeat, min(getattr(gdat, 'fittminm' + namefeat), getattr(gdat, 'trueminm' + namefeat)))
                        setattr(gdat, 'maxm' + namefeat, max(getattr(gdat, 'fittmaxm' + namefeat), getattr(gdat, 'truemaxm' + namefeat)))
                    except:
                        try:
                            setattr(gdat, 'minm' + namefeat, getattr(gdat, 'fittminm' + namefeat))
                            setattr(gdat, 'maxm' + namefeat, getattr(gdat, 'fittmaxm' + namefeat))
                        except:
                            try:
                                setattr(gdat, 'minm' + namefeat, getattr(gdat, 'trueminm' + namefeat))
                                setattr(gdat, 'maxm' + namefeat, getattr(gdat, 'truemaxm' + namefeat))
                            except:
                                pass
    
    # set plot limits
    # temp -- this should be different for each population
    for strg, valu in deepcopy(gdat).__dict__.iteritems():
        if strg.startswith('minm'):
            try:
                minm = getattr(gdat, 'minm' + strg[4:])
                maxm = getattr(gdat, 'maxm' + strg[4:])
                limt = np.array([minm, maxm])
                setattr(gdat, 'limt' + strg[4:], limt)
            except:
                pass

    gdat.boolhash = False
    for strgmodl in gdat.liststrgmodl:
        elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        numbtrap = getattr(gdat, strgmodl + 'numbtrap')
        if numbtrap > 0:
            for l in indxpopl:
                if elemspatevaltype[l] == 'locl':
                    if strgmodl == 'true' and gdat.truenumbelempopl[l] > 0:
                        gdat.boolhash = True
                    if strgmodl == 'fitt' and gdat.fittmaxmnumbelempopl[l] > 0:
                        gdat.boolhash = True
    
    if gdat.rtagmock is not None:
        if gdat.datatype == 'inpt':
            path = gdat.pathoutprtagmock + 'gdatfinlpost'
            booltemp = True
            try:
                gdatmock = readfile(path)
            except:
                print 'Could not read the gdatfinlpost object of the run on the mock data. Setting gdat.rtagmock is None.'
                booltemp = False
                gdat.rtagmock = None

            if booltemp:
                gdat.truenumbtrap = gdatmock.truenumbtrap
                if gdatmock.trueindxpopl != gdat.fittindxpopl:
                    raise Exception('')
                for l in gdat.fittindxpopl:
                    for strgfeat in gdat.fittliststrgfeat[l]:
                        if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':
                            continue

                        if strgfeat[-4:] in gdat.listnamerefr:
                            continue
                        ptfn = getattr(gdatmock, 'trueptfn' + strgfeat + 'pop%d' % l)
                        setattr(gdat, 'trueptfn' + strgfeat + 'pop%d' % l, ptfn)
                gdat.trueliststrgfeat = gdatmock.trueliststrgfeat
    
    if gdat.boolhash:
        gdat.numbprox = 3
        gdat.indxprox = np.arange(gdat.numbprox)
        minmfeatampl = getattr(gdat, 'fittminm' + gdat.fittnamefeatampl[0])
        maxmfeatampl = getattr(gdat, 'fittmaxm' + gdat.fittnamefeatampl[0])
        gdat.binsprox = np.logspace(np.log10(minmfeatampl), np.log10(maxmfeatampl), gdat.numbprox + 1)
        
        # determine the maximum angle at which the contribution of the element will be computed
        if gdat.maxmangleval is None:
            if gdat.exprtype == 'chan':
                gdat.maxmangleval = np.maximum(gdat.maxmangleval, np.array([5., 6., 9.]) / gdat.anglfact)
            elif gdat.exprtype == 'ferm':
                gdat.maxmangleval = np.maximum(gdat.maxmangleval, np.array([7., 9., 15.]) / gdat.anglfact)
            else:
                gdat.maxmangleval = np.empty(gdat.numbprox)
                if gdat.numbpixlfull != 1:
                    for h in gdat.indxprox:
                        if gdat.specfraceval == 0:
                            gdat.maxmangleval[h] = 3. * gdat.maxmgang
                        else:  
                            frac = min(1e-2, gdat.specfraceval * gdat.binsprox[0] / gdat.binsprox[h+1])
                            psfnwdth = retr_psfnwdth(gdat, gdat.psfnexpr, frac)
                            gdat.indxmaxmangl = unravel_index(np.argmax(psfnwdth), psfnwdth.shape)
                            gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
        
        if gdat.commboolelempsfnanyy and gdat.maxmangl - np.amax(gdat.maxmangleval) < 1.1 * np.sqrt(2) * (gdat.maxmlgal - gdat.maxmgangdata):
            print 'gdat.maxmangl'
            print gdat.maxmangl * gdat.anglfact
            print 'gdat.maxmangleval'
            print gdat.maxmangleval * gdat.anglfact
            print 'gdat.maxmlgal'
            print gdat.maxmlgal * gdat.anglfact
            print 'gdat.maxmgangdata'
            print gdat.maxmgangdata * gdat.anglfact
            raise Exception('Angular axis is too short.')

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathprox + 'indxprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * np.amin(gdat.maxmangleval), \
                                                                                                            1e2 * np.amax(gdat.maxmangleval), gdat.numbprox)
    
        if gdat.verbtype > 0 and boolinitsetp:
            print 'Element evaluation will be performed up to'
            if gdat.numbpixlfull != 1:
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
            gdat.indxpixlprox = [[] for h in range(gdat.numbprox)]
            cntrsave = -1.
            # temp
            for j in gdat.indxpixl:
                dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
                dist[j] = 0.
                for h in range(gdat.numbprox):
                    indxpixlproxtemp = np.where(dist < gdat.maxmangleval[h])[0]
                    if indxpixlproxtemp.size > 2e4:
                        indxpixlproxtemp = -1
                        if gdat.maxmangl < np.sqrt(2.) * gdat.maxmgangdata:
                            raise Exception('Angular axis used to interpolate the PSF should be longer.')
                    
                    if indxpixlproxtemp.size < 10:
                        print 'gdat.expo'
                        summgene(gdat.expo)
                        print 'dist'
                        summgene(dist * gdat.anglfact)
                        print 'indxpixlproxtemp.size'
                        print indxpixlproxtemp.size
                        print 'gdat.indxpixl'
                        summgene(gdat.indxpixl)
                        print 'gdat.numbpixl'
                        print gdat.numbpixl
                        print 'gdat.maxmangleval'
                        print gdat.maxmangleval * gdat.anglfact
                        raise Exception('Pixel hash list should not have fewer than 10 pixels.')

                    gdat.indxpixlprox[h].append(indxpixlproxtemp)
                cntrsave = tdpy.util.show_prog(j, gdat.numbpixl, cntrsave)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
        
        gdat.numbpixlprox = np.zeros(gdat.numbprox) 
        for h in range(gdat.numbprox):
            for m in range(len(gdat.indxpixlprox[h])):
                if isinstance(gdat.indxpixlprox[h][m], int):
                    gdat.numbpixlprox[h] += gdat.numbpixl
                else:
                    gdat.numbpixlprox[h] += len(gdat.indxpixlprox[h][m])
            gdat.numbpixlprox[h] /= len(gdat.indxpixlprox[h])
            setattr(gdat, 'numbpixlprox%04d' % h, gdat.numbpixlprox[h])
        if gdat.verbtype > 0:
            print 'gdat.numbpixlprox'
            print gdat.numbpixlprox
        if (gdat.numbpixlprox - np.mean(gdat.numbpixlprox) == 0.).all():
            print 'gdat.lgalgrid'
            summgene(gdat.lgalgrid * gdat.anglfact)
            print 'gdat.bgalgrid'
            summgene(gdat.bgalgrid * gdat.anglfact)
            print 'gdat.maxmgangdata'
            print gdat.maxmgangdata * gdat.anglfact
            if gdat.pixltype == 'cart':
                print 'gdat.numbsidecart'
                print gdat.numbsidecart
            raise Exception('Number of pixels in the hash lists should be different.')

        if gdat.booldiagmode:
            diff = gdat.numbpixlprox - np.roll(gdat.numbpixlprox, 1)
            if not (diff[1:] >= 0).all():
                print 'gdat.numbpixlprox'
                print gdat.numbpixlprox
                print 'np.roll(gdat.numbpixlprox, 1)'
                print np.roll(gdat.numbpixlprox, 1)
                print 'diff'
                print diff
                raise Exception('Number of pixels in the pixel look-up table is wrong.')
    

def setpfinl(gdat, boolinitsetp=False):
    
    gdat.listnamepdir = ['forw', 'reve']
    gdat.listlablpdir = ['f', 'r']
    
    # terms in the log-acceptance probability
    gdat.listnametermlacp = []
    gdat.listlabltermlacp = []
    for l in gdat.fittindxpopl:
        if gdat.fittnumbpopl > 1:
            strgpopl = '%d,' % l
        else:
            strgpopl = ''
        for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
            labl = getattr(gdat, 'labl' + strgcomp)
            gdat.listlabltermlacp += ['$u_{%s%s}$' % (strgpopl, labl)]
    gdat.listnametermlacp += ['lrpp']
    gdat.listlabltermlacp += [u'$\ln P(q)$']
    gdat.listnametermlacp += ['ljcb']
    gdat.listlabltermlacp += [r'$\ln \alpha_j$']
    
    gdat.numbtermlacp = len(gdat.listnametermlacp)
    gdat.indxtermlacp = np.arange(gdat.numbtermlacp)

    if gdat.datatype == 'mock':
        for name in gdat.truenamefixp:
            setattr(gdat, 'true' + name, gdat.truesamp[getattr(gdat, 'trueindxfixp' + name)])

    if gdat.datatype == 'mock' and gdat.truenumbelemtotl != 0:
        for l in gdat.trueindxpopl:
            setattr(gdat, 'truemeanelempop%d' % l, None)
        if gdat.fittlensmodltype != 'none':
            for namecalc in ['delt', 'intg']:
                for nametemp in ['', 'bein']:
                    setattr(gdat, 'truefracsubh%s%s' % (namecalc, nametemp), None)
                    setattr(gdat, 'truemasssubh%s%s' % (namecalc, nametemp), None)
                    setattr(gdat, 'scalfracsubh%s%s' % (namecalc, nametemp), 'self')
                    setattr(gdat, 'scalmasssubh%s%s' % (namecalc, nametemp), 'self')

    # for each parameter in the fitting model, determine if there is a corresponding parameter in the generative model
    for k in gdat.indxvarbscal:
        try:
            temp = getattr(gdat, 'true' + gdat.listnamevarbscal[k])
        except:
            temp = None
        setattr(gdat, 'corr' + gdat.listnamevarbscal[k], temp)
    
    gdat.fittcorrfixp = np.empty(gdat.fittnumbfixp)
    for k in gdat.fittindxfixp:
        try:
            gdat.fittcorrfixp[k] = getattr(gdat, 'true' + gdat.fittnamefixp[k])
        except:
            gdat.fittcorrfixp[k] = None

    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        lpdfprio = [None for l in indxpopl]
        lpdfprioobjt = [None for l in indxpopl]
        lpdfpriointp = [None for l in indxpopl]
        for l in indxpopl:
            liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
            liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
            for strgfeat, strgpdfn in zip(liststrgfeatprio, liststrgpdfnprio):
                if strgpdfn == 'tmplgrad':
                    pdfnpriotemp = np.empty((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                    for k in range(gdat.numbsidecart + 1):
                        for n in range(gdat.numbsidecart + 1):
                            # temp
                            temp = retr_rele(gdat, gdat.truecntp['lens'][0, :, 0], gdat.binslgal[k], gdat.binsbgal[n], 1., \
                                                                                    gdat.trueascadistmeanpop0, gdat.trueacutdistmeanpop0, gdat.indxpixl)
    
                            #temp /= np.amax(temp)
                            pdfnpriotemp[k, n] = temp**gdat.relnpowr
                    lpdfprio, lpdfprioobjt = retr_spatprio(gdat, pdfnpriotemp)
                    lpdfpriointp = lpdfprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
        setattr(gdat, strgmodl + 'lpdfprio', lpdfprio)
        setattr(gdat, strgmodl + 'lpdfprioobjt', lpdfprioobjt)
        setattr(gdat, strgmodl + 'lpdfpriointp', lpdfpriointp)
        
    # find the pixels at which data count maps have local maxima
    if gdat.pixltype == 'cart':
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                # temp
                gdat.indxxdatmaxm, gdat.indxydatmaxm = tdpy.util.retr_indximagmaxm(gdat.cntpdatacart[i, :, m])

    # sanity checks
    # temp
    if (np.fabs(gdat.cntpdata - np.round(gdat.cntpdata)) > 1e-3).any() and boolinitsetp:
        print 'Fractional counts!'

    if np.amin(gdat.cntpdata) < 0. and boolinitsetp:
        print 'Negative counts!'


def retr_gradmaps(gdat, maps):
    
    # temp -- this does not work with vanishing exposure
    maps = maps.reshape((gdat.numbsidecart, gdat.numbsidecart))
    grad = np.dstack((np.gradient(maps, gdat.sizepixl, axis=0), np.gradient(maps, gdat.sizepixl, axis=1))).reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    grad = grad.reshape((gdat.numbpixlcart, 2))

    return grad


def retr_spatmean(gdat, inpt, boolcntp=False):
    
    listspatmean = [[] for b in gdat.indxspatmean]
    listspatstdv = [[] for b in gdat.indxspatmean]
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        if boolcntp:
            cntp = inpt[gdat.listindxcubespatmean[b]]
        else:
            cntp = inpt[gdat.listindxcubespatmean[b]] * gdat.expo[gdat.listindxcubespatmean[b]] * gdat.apix
            if gdat.enerdiff:
                cntp *= gdat.deltener[:, None, None]
        spatmean = np.mean(np.sum(cntp, 2), axis=1) / gdat.expototlmean / gdat.apix
        spatstdv = np.sqrt(np.sum(cntp, axis=(1, 2))) / gdat.numbdata / gdat.expototlmean / gdat.apix
        if gdat.enerdiff:
            spatmean /= gdat.deltener
            spatstdv /= gdat.deltener
        listspatmean[b] = spatmean
        listspatstdv[b] = spatstdv

    return listspatmean, listspatstdv


def retr_rele(gdat, maps, lgal, bgal, defs, asca, acut, indxpixlelem, absv=True, cntpmodl=None):
    
    grad = retr_gradmaps(gdat, maps)
        
    defl = retr_defl(gdat, indxpixlelem, lgal, bgal, defs, asca=asca, acut=acut)

    prod = grad * defl
    if cntpmodl is not None:
        prod /= cntpmodl[:, None]
    dotstemp = np.sum(prod, 1)
    if absv:
        dotstemp = np.fabs(dotstemp)
    else:
        dotstemp = dotstemp
    
    dots = np.mean(dotstemp)
    
    return dots


def timefunc(name, func, *args, **kwargs):

    func(*args, **kwargs)
    minm = 1e10
    for k in range(5):
        timeinit = time.clock()
        func(*args, **kwargs)
        timefinl = time.clock()
        timediff = timefinl - timeinit
        minm = min(timediff, minm)
    
    print 'Timing %s...' % name
    print '%.3g ms' % (1e3 * (timefinl - timeinit))


def retr_listconc(thislist):
    
    thislistconc = []
    for listsubb in thislist:
        for strg in listsubb:
            if not strg in thislistconc:
                thislistconc.append(strg)
    
    return thislistconc


def retr_datatick(gdat):

    # data count limits
    gdat.maxmcntpdata = 0
    gdat.maxmcntpdata = max(0.5, np.amax(np.sum(gdat.cntpdata, 2)))
    gdat.maxmcntpdata = np.amax(gdat.maxmcntpdata)
    gdat.maxmcntpmodl = gdat.maxmcntpdata
    gdat.minmcntpdata = max(0.5, 1e-4 * gdat.maxmcntpdata)
    gdat.minmcntpmodl = 1e-1 * gdat.minmcntpdata
    gdat.maxmcntpresi = np.ceil(gdat.maxmcntpdata * 0.1)
    gdat.minmcntpresi = -gdat.maxmcntpresi

    # tick labels for count maps
    liststrgcbar = ['cntpdata', 'cntpresi', 'cntpmodl']
    for strgcbar in liststrgcbar:
        retr_ticklabl(gdat, strgcbar)


def retr_ticklabl(gdat, strgcbar):

    scal = getattr(gdat, 'scal' + strgcbar)
    
    minm = getattr(gdat, 'minm' + strgcbar)
    maxm = getattr(gdat, 'maxm' + strgcbar)
    retr_axis(gdat, strgcbar, minm, maxm, gdat.numbtickcbar - 1, scal=scal)

    minmscal = minm
    if scal == 'asnh':
        minmscal = np.arcsinh(minmscal)
    if scal == 'logt':
        minmscal = np.log10(minmscal)
    maxmscal = maxm
    if scal == 'asnh':
        maxmscal = np.arcsinh(maxmscal)
    if scal == 'logt':
        maxmscal = np.log10(maxmscal)

    tickscal = np.linspace(minmscal, maxmscal, gdat.numbtickcbar)
    labl = np.empty(gdat.numbtickcbar, dtype=object)
    tick = np.copy(tickscal)
    for k in range(gdat.numbtickcbar):
        if scal == 'asnh':
            tick[k] = np.sinh(tickscal[k])
        elif scal == 'logt':
            tick[k] = 10**(tickscal[k])

        # avoid very small, but nonzero central values in the residual count color maps
        if strgcbar == 'cntpresi' and np.fabs(tick[k]) < 1e-5:
            tick[k] = 0.

        if strgcbar == 'cntpdata' and np.amax(tick) > 1e3:
            labl[k] = '%d' % tick[k]
        else:
            labl[k] = '%.3g' % tick[k]
    
    setattr(gdat, 'tick' + strgcbar, tickscal)
    setattr(gdat, 'minmscal' + strgcbar, minmscal)
    setattr(gdat, 'maxmscal' + strgcbar, maxmscal)
    setattr(gdat, 'labl' + strgcbar, labl)
    

def retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn, strgmome='pmea', indxvarb=None, indxlist=None):
    
    if strgvarb.startswith('cntpdata'):
        varb = getattr(gdat, strgvarb)
    elif strgvarb.startswith('histcntpdata'):
        varb = getattr(gdat, strgvarb)
    else:
        if strgmodl == 'true':
            varb = getattr(gdat, 'true' + strgvarb)
        if strgmodl == 'fitt':
            if strgstat == 'this':
                if strgmome == 'errr':
                    varb = getattr(gdatmodi, strgstat + 'errr' + strgvarb)
                else:
                    varb = getattr(gdatmodi, strgstat + strgvarb)
            if strgstat == 'pdfn':
                varb = getattr(gdat, strgmome + strgpdfn + strgvarb)

    if indxlist is not None:
        varb = varb[indxlist]

    if indxvarb is not None:
        if strgmome == 'errr':
            varb = varb[[slice(None)] + indxvarb]
        else:
            varb = varb[indxvarb]

    return np.copy(varb)


def retr_indxsamp(gdat, strgmodl='fitt', init=False):
    
    dicttemp = {}
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    if lensmodltype != 'none' or hostemistype != 'none':
        numbsersfgrd = getattr(gdat, strgmodl + 'numbsersfgrd')

    if init:
        # transdimensional element populations
        backtype = getattr(gdat, strgmodl + 'backtype')
        
        numbback = 0
        indxback = []
        for c in range(len(backtype)):
            if isinstance(backtype[c], str):
                if backtype[c].startswith('bfunfour') or backtype[c].startswith('bfunwfou'):
                    namebfun = backtype[c][:8]
                    ordrexpa = int(backtype[c][8:])
                    numbexpa = 4 * ordrexpa**2
                    indxexpa = np.arange(numbexpa)
                    del backtype[c]
                    for k in indxexpa:
                        backtype.insert(c+k, namebfun + '%04d' % k)
        numbback = len(backtype)
        indxback = np.arange(numbback)
        numbbacktotl = np.sum(numbback)
        indxbacktotl = np.arange(numbbacktotl)
        numbpopl = len(elemtype)
        
        indxpopl = np.arange(numbpopl)
        if lensmodltype != 'none' or hostemistype != 'none':
            indxsersfgrd = []
            indxsersfgrd = np.arange(numbsersfgrd)

        # feature used to model the amplitude of elements
        namefeatampl = [[] for l in indxpopl]
        indxcompampl = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l] == 'lghtpntspuls':
                namefeatampl[l] = 'per0'
                indxcompampl[l] = 2
            elif elemtype[l] == 'lghtpntsagnntrue':
                namefeatampl[l] = 'lum0'
                indxcompampl[l] = 2
            elif elemtype[l].startswith('lghtline'):
                namefeatampl[l] = 'flux'
                indxcompampl[l] = 1
            elif elemtype[l].startswith('lghtpnts'):
                namefeatampl[l] = 'flux'
                indxcompampl[l] = 2
            elif elemtype[l].startswith('lghtgausbgrd'):
                namefeatampl[l] = 'flux'
                indxcompampl[l] = 2
            if elemtype[l] == 'lens':
                namefeatampl[l] = 'defs'
                indxcompampl[l] = 2
            if elemtype[l].startswith('clus'):
                namefeatampl[l] = 'nobj'
                indxcompampl[l] = 2
            if len(namefeatampl[l]) == 0:
                raise Exception('Amplitude feature undefined.')
    else:
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
        namefeatampl = getattr(gdat, strgmodl + 'namefeatampl') 
        indxcompampl = getattr(gdat, strgmodl + 'indxcompampl') 
        maxmnumbelem = getattr(gdat, strgmodl + 'maxmnumbelem') 
        spatdisttype = getattr(gdat, strgmodl + 'spatdisttype') 
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        indxback = getattr(gdat, strgmodl + 'indxback')
        listnameback = getattr(gdat, strgmodl + 'listnameback') 
        
        if strgmodl == 'true':
            numbelempopl = np.zeros(numbpopl)
            for l in indxpopl:
                numbelempopl[l] += getattr(gdat, 'truenumbelempop%d' % l)
         
        numbelemzero = [[] for l in indxpopl]
        for l in indxpopl:
            numbelemzero[l] = 0

        # element setup
        ## flag to calculate the kernel approximation errors
        calcerrr = [[] for l in indxpopl]
        for l in indxpopl:
            if elemspatevaltype[l] != 'full' and gdat.numbpixlfull < 1e5:
                # temp
                calcerrr[l] = False
            else:
                calcerrr[l] = False
        
        # total np.maximum number of elements
        maxmnumbelempopl = np.zeros(numbpopl, dtype=int)
        for l in indxpopl:
            maxmnumbelempopl[l] += np.sum(maxmnumbelem[l])
        maxmnumbelemtotl = np.sum(maxmnumbelempopl) 

        ## sorting feature
        namefeatsort = [[] for l in indxpopl]
        for l in indxpopl:
            # feature to be used to sort elements
            if elemtype[l].startswith('lght'):
                namefeatsort[l] = 'flux'
            if elemtype[l] == 'lens':
                namefeatsort[l] = 'defs'
            if elemtype[l].startswith('clus'):
                namefeatsort[l] = 'nobj'
    
        ## selection feature
        listnamefeatsele = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght'):
                listnamefeatsele[l] = ['flux']
            if elemtype[l] == 'lens':
                listnamefeatsele[l] = ['defs', 'mcut', 'rele']
            if elemtype[l].startswith('clus'):
                listnamefeatsele[l] = ['nobj']
    
        ## label extensions
        lablelemextn = [[] for l in indxpopl]
        for l in indxpopl:
            if gdat.numbgrid > 1:
                if elemtype[l] == 'lghtpnts':
                    lablelemextn[l] = r'\rm{fps}'
                if elemtype[l] == 'lghtgausbgrd':
                    lablelemextn[l] = r'\rm{bgs}'
            else:
                if elemtype[l].startswith('lghtpntspuls'):
                    lablelemextn[l] = r'\rm{pul}'
                if elemtype[l].startswith('lghtpntsagnn'):
                    lablelemextn[l] = r'\rm{agn}'
                elif elemtype[l] == 'lghtpnts':
                    lablelemextn[l] = r'\rm{pts}'
            if elemtype[l] == 'lens':
                lablelemextn[l] = r'\rm{sub}'
            if elemtype[l].startswith('clus'):
                lablelemextn[l] = r'\rm{cls}'
            if elemtype[l].startswith('lghtline'):
                lablelemextn[l] = r'\rm{lin}'
    
        indxpoplgrid = [[] for y in gdat.indxgrid]
        for y in gdat.indxgrid: 
            for indx, elemtypetemp in enumerate(elemtype):
                # foreground grid (image plane) -- the one np.where the data is measured
                if y == 0:
                    if elemtypetemp.startswith('lght') and not elemtypetemp.endswith('bgrd') or elemtypetemp.startswith('clus'):
                        indxpoplgrid[y].append(indx)
                # foreground mass grid
                if y == 1:
                    if elemtypetemp.startswith('lens'):
                        indxpoplgrid[y].append(indx)
                # background grid (source plane)
                if y == 2:
                    if elemtypetemp.endswith('bgrd'):
                        indxpoplgrid[y].append(indx)
        
        indxgridpopl = [[] for l in indxpopl]
        for l in indxpopl:
            for y in gdat.indxgrid:
                if l in indxpoplgrid[y]:
                    indxgridpopl[l] = y
    
        calcelemsbrt = False
        for l in indxpopl:
            if elemtype[l].startswith('lghtpnts'):
                calcelemsbrt = True
    
        if 'lghtgausbgrd' in elemtype:
            calcelemsbrtbgrd = True
        else:
            calcelemsbrtbgrd = False

        if 'lens' in elemtype:
            calcelemdefl = True
        else:
            calcelemdefl = False

        ## element Boolean flags
        boolelemlght = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght'):
                boolelemlght[l] = True
            else:
                boolelemlght[l] = False
        boolelemlghtanyy = True in boolelemlght

        boolelemspat = [[] for l in indxpopl]
        for l in indxpopl:
            if not elemtype[l].startswith('lghtline'):
                boolelemspat[l] = True
            else:
                boolelemspat[l] = False
        boolelemspatanyy = True in boolelemspat

        boolelemlghtspat = [[] for l in indxpopl]
        for l in indxpopl:
            
            if boolelemspat[l] and boolelemlght[l]:
                boolelemlghtspat[l] = True
            else:
                boolelemlghtspat[l] = False

        boolelemlens = False
        for l in indxpopl:
            if elemtype[l].startswith('lens'):
                boolelemlens = True
        
        boolelemsbrtdfnc = [[] for l in indxpopl]
        for l in indxpopl:
            if maxmnumbelempopl[l] > 0 and (elemtype[l].startswith('lght') and not elemtype[l].endswith('bgrd') or elemtype[l].startswith('clus')):
                boolelemsbrtdfnc[l] = True
            else:
                boolelemsbrtdfnc[l] = False
        boolelemsbrtdfncanyy = True in boolelemsbrtdfnc

        boolelemdeflsubh = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l] == 'lens':
                boolelemdeflsubh[l] = True
            else:
                boolelemdeflsubh[l] = False
        boolelemdeflsubhanyy = True in boolelemdeflsubh

        boolelemsbrtextsbgrd = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght') and elemtype[l].endswith('bgrd'):
                boolelemsbrtextsbgrd[l] = True
            else:
                boolelemsbrtextsbgrd[l] = False
        boolelemsbrtextsbgrdanyy = True in boolelemsbrtextsbgrd
        
        if boolelemsbrtextsbgrdanyy:
            indxpopllens = 1
        else:
            indxpopllens = 0

        boolelemsbrtpnts = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght') and elemtype[l] != 'lghtline' or elemtype[l] == 'clus':
                boolelemsbrtpnts[l] = True
            else:
                boolelemsbrtpnts[l] = False
        boolelemsbrtpntsanyy = True in boolelemsbrtpnts

        # temp -- because there is currently no extended source
        boolelemsbrt = boolelemsbrtdfnc
    
        boolelempsfn = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lghtpnts') or elemtype[l] == 'clus':
                boolelempsfn[l] = True
            else:
                boolelempsfn[l] = False
        boolelempsfnanyy = True in boolelempsfn
        
        spectype = [[] for l in indxpopl]
        for l in indxpopl:
            if boolelemlght[l]:
                spectype[l] = 'powr'
            else:
                spectype[l] = 'none'
        setp_varbvalu(gdat, 'spectype', spectype, strgmodl=strgmodl)
    
        minmgwdt = 2. * gdat.sizepixl
        maxmgwdt = gdat.maxmgangdata / 4.
        setp_varblimt(gdat, 'gwdt', [minmgwdt, maxmgwdt], strgmodl=strgmodl)
    
        if boolelemlghtanyy:
            # flux
            if gdat.exprtype == 'ferm':
                minmflux = 1e-9
                maxmflux = 1e-6
            if gdat.exprtype == 'tess':
                minmflux = 1.
                maxmflux = 1e3
            if gdat.exprtype == 'chan':
                if gdat.anlytype == 'spec':
                    minmflux = 1e4
                    maxmflux = 1e7
                else:
                    minmflux = 3e-9
                    maxmflux = 1e-6
            if gdat.exprtype == 'sdyn':
                minmflux = 0.1
                maxmflux = 100.
            if gdat.exprtype == 'hubb':
                minmflux = 1e-20
                maxmflux = 1e-17
            if gdat.exprtype == 'fire':
                minmflux = 1e-20
                maxmflux = 1e-17
            setp_varblimt(gdat, 'flux', [minmflux, maxmflux], strgmodl=strgmodl)
            
            if gdat.exprtype == 'ferm':
                setp_varblimt(gdat, 'fluxdistbrek', [3e-9, 1e-6], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'fluxdistsloplowr', [0.5, 3.], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'fluxdistslopuppr', [0.5, 3.], popl=l, strgmodl=strgmodl)
            
            if gdat.enerbins:
                ### spectral parameters
                if gdat.exprtype == 'ferm':
                    sind = [1., 3.]
                    minmsind = 1.
                    maxmsind = 3.
                if gdat.exprtype == 'chan':
                    minmsind = 0.4
                    maxmsind = 2.4
                    sind = [0.4, 2.4]
                if gdat.exprtype == 'hubb':
                    minmsind = 0.5
                    maxmsind = 2.5
                    sind = [0.4, 2.4]
                if gdat.exprtype != 'fire':
                    setp_varblimt(gdat, 'sind', [minmsind, maxmsind], strgmodl=strgmodl)
                    setp_varblimt(gdat, 'curv', [-1., 1.], strgmodl=strgmodl)
                    setp_varblimt(gdat, 'expc', [0.1, 10.], strgmodl=strgmodl)
                    setp_varblimt(gdat, 'sinddistmean', sind, popl='full', strgmodl=strgmodl)
                    #### standard deviations should not be too small
                    setp_varblimt(gdat, 'sinddiststdv', [0.3, 2.], popl='full', strgmodl=strgmodl)
                    setp_varblimt(gdat, 'curvdistmean', [-1., 1.], popl='full', strgmodl=strgmodl)
                    setp_varblimt(gdat, 'curvdiststdv', [0.1, 1.], popl='full', strgmodl=strgmodl)
                    setp_varblimt(gdat, 'expcdistmean', [1., 8.], popl='full', strgmodl=strgmodl)
                    setp_varblimt(gdat, 'expcdiststdv', [0.01 * gdat.maxmener, gdat.maxmener], popl='full', strgmodl=strgmodl)
                    for i in gdat.indxenerinde:
                        setp_varblimt(gdat, 'sindcolr0001', [-2., 6.], strgmodl=strgmodl)
                        setp_varblimt(gdat, 'sindcolr0002', [0., 8.], strgmodl=strgmodl)
                        #setp_varblimt(gdat, 'sindcolr%04d' % i, [-5., 10.], strgmodl=strgmodl)
        
        for l in indxpopl:
            if elemtype[l] == 'lghtpntspuls':
                setp_varblimt(gdat, 'gang', [1e-1 * gdat.sizepixl, gdat.maxmgangdata], strgmodl=strgmodl)
                setp_varblimt(gdat, 'geff', [0., 0.4], strgmodl=strgmodl)
                setp_varblimt(gdat, 'dglc', [10., 3e3], strgmodl=strgmodl)
                setp_varblimt(gdat, 'phii', [0., 2. * np.pi], strgmodl=strgmodl)
                setp_varblimt(gdat, 'thet', [0., np.pi], strgmodl=strgmodl)
                setp_varblimt(gdat, 'per0distmean', [5e-4, 1e1], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'magfdistmean', [1e7, 1e16], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'per0diststdv', [1e-2, 1.], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'magfdiststdv', [1e-2, 1.], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'gangdistslop', [0.5, 4.], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'dglcdistslop', [0.5, 2.], popl=l, strgmodl=strgmodl)
                # temp -- minima and maxima should be population dependent...
    
                setp_varblimt(gdat, 'spatdistcons', [1e-4, 1e-2], popl='full')
                setp_varblimt(gdat, 'bgaldistscal', [0.5 / gdat.anglfact, 5. / gdat.anglfact], popl='full', strgmodl=strgmodl)
            if elemtype[l] == 'lghtpntsagnntrue':
                setp_varblimt(gdat, 'dlos', [1e7, 1e9], strgmodl=strgmodl)
                setp_varblimt(gdat, 'dlosdistslop', [-0.5, -3.], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'lum0', [1e43, 1e46], strgmodl=strgmodl)
                setp_varblimt(gdat, 'lum0distbrek', [1e42, 1e46], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'lum0distsloplowr', [0.5, 3.], popl=l, strgmodl=strgmodl)
                setp_varblimt(gdat, 'lum0distslopuppr', [0.5, 3.], popl=l, strgmodl=strgmodl)
        
        # construct background surface brightness templates from the user input
        lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
        if lensmodltype != 'none' or hostemistype != 'none':
            indxsersfgrd = getattr(gdat, strgmodl + 'indxsersfgrd')
        backtype = getattr(gdat, strgmodl + 'backtype')
        numbback = getattr(gdat, strgmodl + 'numbback')
        sbrtbacknorm = [[] for c in indxback]
        unifback = np.ones(numbback, dtype=bool)
        for c in indxback:
            sbrtbacknorm[c] = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            if backtype[c] == 'data':
                sbrtbacknorm[c] = np.copy(gdat.sbrtdata)
                sbrtbacknorm[c][np.where(sbrtbacknorm[c] == 0.)] = 1e-100
            elif isinstance(backtype[c], float):
                sbrtbacknorm[c] = np.zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + backtype[c]
            elif isinstance(backtype[c], list) and isinstance(backtype[c], float):
                sbrtbacknorm[c] = retr_spec(gdat, np.array([backtype[c]]), sind=np.array([backtype[c]]))[:, 0, None, None]
            elif isinstance(backtype[c], np.ndarray) and backtype[c].ndim == 1:
                sbrtbacknorm[c] = np.zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + backtype[c][:, None, None]
            elif backtype[c].startswith('bfunfour') or backtype[c].startswith('bfunwfou'):
                namebfun = getattr(gdat, strgmodl + 'namebfun')
                ordrexpa = getattr(gdat, strgmodl + 'ordrexpa')
                numbexpa = getattr(gdat, strgmodl + 'numbexpa')
                indxexpatemp = int(backtype[c][8:]) 
                indxterm = indxexpatemp // ordrexpa**2
                indxexpaxdat = (indxexpatemp % ordrexpa**2) // ordrexpa + 1
                indxexpaydat = (indxexpatemp % ordrexpa**2) % ordrexpa + 1
                if namebfun == 'bfunfour':
                    ampl = 1.
                    func = gdat.meanbgalcart 
                if namebfun == 'bfunwfou':
                    functemp = np.exp(-0.5 * (gdat.meanbgalcart / (1. / gdat.anglfact))**2)
                    ampl = np.sqrt(functemp)
                    func = functemp
                argslgal = 2. * np.pi * indxexpaxdat * gdat.meanlgalcart / gdat.maxmgangdata
                argsbgal = 2. * np.pi * indxexpaydat * func / gdat.maxmgangdata
                if indxterm == 0:
                    termfrst = np.sin(argslgal)
                    termseco = ampl * np.sin(argsbgal)
                if indxterm == 1:
                    termfrst = np.sin(argslgal)
                    termseco = ampl * np.cos(argsbgal)
                if indxterm == 2:
                    termfrst = np.cos(argslgal)
                    termseco = ampl * np.sin(argsbgal)
                if indxterm == 3:
                    termfrst = np.cos(argslgal)
                    termseco = ampl * np.cos(argsbgal)
                sbrtbacknorm[c] = (termfrst[None, :] * termseco[:, None]).flatten()[None, :, None] * np.ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                
            else:
                path = gdat.pathinpt + backtype[c]
                sbrtbacknorm[c] = astropy.io.fits.getdata(path)
                
                if gdat.pixltype == 'cart':
                    if not gdat.forccart:
                        if sbrtbacknorm[c].shape[2] != gdat.numbsidecart:
                            print 'gdat.numbsidecart'
                            print gdat.numbsidecart
                            print 'sbrtbacknorm[c]'
                            summgene(sbrtbacknorm[c])
                            raise Exception('Provided background template must have the chosen image dimensions.')
                    
                    sbrtbacknorm[c] = sbrtbacknorm[c].reshape((sbrtbacknorm[c].shape[0], -1, sbrtbacknorm[c].shape[-1]))
        
                if gdat.pixltype == 'cart' and gdat.forccart:
                    sbrtbacknormtemp = np.empty((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                    for i in gdat.indxenerfull:
                        for m in gdat.indxevttfull:
                            sbrtbacknormtemp[i, :, m] = tdpy.util.retr_cart(sbrtbacknorm[c][i, :, m], \
                                                    numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                    minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                    minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).flatten()
                    sbrtbacknorm[c] = sbrtbacknormtemp

            # determine spatially uniform background templates
            for i in gdat.indxenerfull:
                for m in gdat.indxevttfull:
                    if np.std(sbrtbacknorm[c][i, :, m]) > 1e-6:
                        unifback[c] = False

        boolzero = True
        boolbfun = False
        for c in indxback:
            if np.amin(sbrtbacknorm[c]) < 0. and isinstance(backtype[c], str) and not backtype[c].startswith('bfun'):
                booltemp = False
                raise Exception('Background templates must be positive-definite everynp.where.')
        
            if not np.isfinite(sbrtbacknorm[c]).all():
                raise Exception('Background template is not finite.')

            if np.amin(sbrtbacknorm[c]) > 0. or backtype[c] == 'data':
                boolzero = False
            
            if isinstance(backtype[c], str) and backtype[c].startswith('bfun'):
                boolbfun = True
        
        if boolzero and not boolbfun:
            raise Exception('At least one background template must be positive everynp.where.')
       
        if maxmnumbelemtotl > 0 and boolelempsfnanyy:
            if hostemistype != 'none' or not unifback.all():
                psfnevaltype = 'full'
            else:
                psfnevaltype = 'kern'
        else:
            if hostemistype != 'none' or not unifback.all():
                psfnevaltype = 'conv'
            else:
                psfnevaltype = 'none'
        
        setp_varbvalu(gdat, 'psfnevaltype', psfnevaltype, strgmodl=strgmodl)
        psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')

        ### PSF model
        if psfnevaltype != 'none':
            psfntype = getattr(gdat, strgmodl + 'psfntype')
            
            print 'psfntype'
            print psfntype
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
            
            numbpsfptotl = numbpsfpform
            
            if gdat.psfninfoprio:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        meansigc = gdat.psfpexpr[i * numbpsfptotl + m * numbpsfptotl * gdat.numbener]
                        stdvsigc = meansigc * 0.1
                        setp_varblimt(gdat, 'sigcen%02devt%d' % (i, m), [meansigc, stdvsigc], typelimt='meanstdv', strgmodl=strgmodl)
                        if psfntype == 'doubking' or psfntype == 'singking':
                            meangamc = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 1]
                            stdvgamc = meangamc * 0.1
                            setp_varblimt(gdat, 'gamcen%02devt%d' % (i, m), [meangamc, stdvgamc], typelimt='meanstdv', strgmodl=strgmodl)
                            if psfntype == 'doubking':
                                meansigt = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 2]
                                stdvsigt = meansigt * 0.1
                                setp_varblimt(gdat, 'sigten%02devt%d' % (i, m), [meansigt, stdvsigt], typelimt='meanstdv', strgmodl=strgmodl)
                                meangamt = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 3]
                                stdvgamt = meangamt * 0.1
                                setp_varblimt(gdat, 'gamten%02devt%d' % (i, m), [meangamt, stdvgamt], typelimt='meanstdv', strgmodl=strgmodl)
                                meanpsff = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 4]
                                stdvpsff = meanpsff * 0.1
                                setp_varblimt(gdat, 'psffen%02devt%d' % (i, m), [meanpsff, stdvpsff], typelimt='meanstdv', strgmodl=strgmodl)
            else:
                if gdat.exprtype == 'sdyn':
                    minmsigm = 0.01 / gdat.anglfact
                    maxmsigm = 0.1 / gdat.anglfact
                if gdat.exprtype == 'ferm':
                    minmsigm = 0.1
                    maxmsigm = 10.
                if gdat.exprtype == 'hubb':
                    minmsigm = 0.01 / gdat.anglfact
                    maxmsigm = 0.1 / gdat.anglfact
                if gdat.exprtype == 'chan':
                    minmsigm = 0.1 / gdat.anglfact
                    maxmsigm = 2. / gdat.anglfact
                minmgamm = 1.5
                maxmgamm = 20.
                setp_varblimt(gdat, 'sigc', [minmsigm, maxmsigm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'sigt', [minmsigm, maxmsigm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'gamc', [minmgamm, maxmgamm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'gamt', [minmgamm, maxmgamm], ener='full', evtt='full', strgmodl=strgmodl)
            setp_varblimt(gdat, 'psff', [0., 1.], ener='full', evtt='full', strgmodl=strgmodl)
 
        specback = getattr(gdat, strgmodl + 'specback')
        
        spectype = getattr(gdat, strgmodl + 'spectype')

        # background
        ## number of background parameters
        numbbacp = 0
        for c in indxback:
            if specback[c]:
                numbbacp += 1
            else:
                numbbacp += gdat.numbener
   
        ## background parameter indices
        indxbackbacp = np.zeros(numbbacp, dtype=int)
        indxenerbacp = np.zeros(numbbacp, dtype=int)
        cntr = 0
        for c in indxback:
            if specback[c]:
                indxbackbacp[cntr] = c
                cntr += 1
            else:
                for i in gdat.indxener:
                    indxenerbacp[cntr] = i
                    indxbackbacp[cntr] = c
                    cntr += 1
        
        indxbacpback = [[] for c in indxback]
        for c in indxback:
            indxbacpback[c] = np.where((indxbackbacp == c))[0]
                
        # features which correlate with significance
        liststrgfeatsign = [[] for l in indxpopl]
        for l in indxpopl:
            for namefeat in listnamefeatsele[l]:
                liststrgfeatsign[l] += [namefeat]
        
        if gdat.verbtype > 0:
            if strgmodl == 'true':
                print 'Building elements for the true model...'
            else:
                print 'Building elements for the fitting model...'

        liststrgcomp = [[] for l in indxpopl]
        listscalcomp = [[] for l in indxpopl]
        liststrgfeatscal = [{} for l in indxpopl]
        for l in indxpopl:
            
            if elemtype[l].startswith('lghtline'):
                liststrgcomp[l] = ['elin']
                listscalcomp[l] = ['logt']
            elif spatdisttype[l] == 'diskscal':
                liststrgcomp[l] = ['lgal', 'bgal']
                listscalcomp[l] = ['self', 'dexpscal']
            elif spatdisttype[l] == 'gangexpo':
                liststrgcomp[l] = ['gang', 'aang']
                listscalcomp[l] = ['expo', 'self']
            elif spatdisttype[l] == 'glc3':
                liststrgcomp[l] = ['dglc', 'thet', 'phii']
                listscalcomp[l] = ['powrslop', 'self', 'self']
            else:
                liststrgcomp[l] = ['lgal', 'bgal']
                listscalcomp[l] = ['self', 'self']
            
            # amplitude
            if elemtype[l] == 'lghtpntsagnntrue':
                liststrgcomp[l] += ['lum0']
                listscalcomp[l] += ['dpowslopbrek']
            elif elemtype[l] == 'lghtpntspuls':
                liststrgcomp[l] += ['per0']
                listscalcomp[l] += ['lnormeanstdv']
            elif elemtype[l].startswith('lght'):
                liststrgcomp[l] += ['flux']
                listscalcomp[l] += [getattr(gdat, strgmodl + 'fluxdisttype')[l]]
            elif elemtype[l] == 'lens':
                liststrgcomp[l] += ['defs']
                listscalcomp[l] += ['powrslop']
            elif elemtype[l].startswith('clus'):
                liststrgcomp[l] += ['nobj']
                listscalcomp[l] += ['powrslop']
           
            # shape
            if elemtype[l] == 'lghtgausbgrd' or elemtype[l] == 'clusvari':
                liststrgcomp[l] += ['gwdt']
                listscalcomp[l] += ['powrslop']
            if elemtype[l] == 'lghtlinevoig':
                liststrgcomp[l] += ['sigm']
                listscalcomp[l] += ['logt']
                liststrgcomp[l] += ['gamm']
                listscalcomp[l] += ['logt']
            
            # others
            if elemtype[l] == 'lghtpntspuls':
                liststrgcomp[l] += ['magf']
                listscalcomp[l] += ['lnormeanstdv']
                liststrgcomp[l] += ['geff']
                listscalcomp[l] += ['self']
            elif elemtype[l] == 'lghtpntsagnntrue':
                liststrgcomp[l] += ['dlos']
                listscalcomp[l] += ['powrslop']

            if gdat.numbener > 1 and elemtype[l].startswith('lghtpnts'):
                if spectype[l] == 'colr':
                    for i in gdat.indxener:
                        if i == 0:
                            continue
                        liststrgcomp[l] += ['sindcolr%04d' % i]
                        listscalcomp[l] += ['self']
                else:
                    liststrgcomp[l] += ['sind']
                    listscalcomp[l] += ['self']
                    if spectype[l] == 'curv':
                        liststrgcomp[l] += ['curv']
                        listscalcomp[l] += ['self']
                    if spectype[l] == 'expc':
                        liststrgcomp[l] += ['expc']
                        listscalcomp[l] += ['self']
            if elemtype[l] == 'lens':
                if gdat.variasca:
                    liststrgcomp[l] += ['asca']
                    listscalcomp[l] += ['self']
                if gdat.variacut:
                    liststrgcomp[l] += ['acut']
                    listscalcomp[l] += ['self']
            
            for scaltype in gdat.listscaltype:
                liststrgfeatscal[l][scaltype] = []
                for k, strgcomp in enumerate(liststrgcomp[l]):
                    if scaltype == listscalcomp[l][k]:
                        liststrgfeatscal[l][scaltype].append(strgcomp)

        # variables for which whose marginal distribution and pair-correlations will be plotted
        liststrgfeatodim = [[] for l in indxpopl]
        for l in indxpopl:
            liststrgfeatodim[l] = deepcopy(liststrgcomp[l])
            liststrgfeatodim[l] += ['deltllik']
            if gdat.numbpixlfull > 1:
                if not 'lgal' in liststrgfeatodim[l]:
                    liststrgfeatodim[l] += ['lgal']
                if not 'bgal' in liststrgfeatodim[l]:
                    liststrgfeatodim[l] += ['bgal']
                if not 'gang' in liststrgfeatodim[l]:
                    liststrgfeatodim[l] += ['gang']
                if not 'aang' in liststrgfeatodim[l]:
                    liststrgfeatodim[l] += ['aang']
            if elemtype[l].startswith('lght'):
                liststrgfeatodim[l] += ['cnts']
                if gdat.exprtype == 'ferm':
                    liststrgfeatodim[l] + ['sbrt0018']
                
            if elemtype[l] == 'lghtpntsagnntrue':
                liststrgfeatodim[l] += ['reds']
                liststrgfeatodim[l] += ['lumi']
                liststrgfeatodim[l] += ['flux']
            if elemtype[l] == 'lghtpntspuls':
                liststrgfeatodim[l] += ['lumi']
                liststrgfeatodim[l] += ['flux']
                liststrgfeatodim[l] += ['mass']
                liststrgfeatodim[l] += ['dlos']
            if elemtype[l] == 'lens':
                liststrgfeatodim[l] += ['mcut', 'diss', 'rele', 'reln', 'relk', 'relf', 'relm', 'reld', 'relc']
        
            if strgmodl == 'fitt':
                for q in gdat.indxrefr: 
                    if gdat.fittnamefeatampl[l] in gdat.refrliststrgfeat[q]:
                        liststrgfeatodim[l].append('aerr' + gdat.listnamerefr[q])

        # add reference element features that are not available in the fitting model
        gdat.refrliststrgfeatonly = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
        if strgmodl == 'fitt':
            gdat.fittliststrgfeatextr = [[] for l in gdat.fittindxpopl]
            for q in gdat.indxrefr: 
                if gdat.refrnumbelempopl[q] == 0:
                    continue
                for name in gdat.refrliststrgfeat[q]:
                    for l in indxpopl:
                        if elemtype[l].startswith('lght') and (name == 'defs' or name == 'acut' or name == 'asca' or name == 'mass'):
                            continue
                        if elemtype[l] == ('lens') and (name == 'cnts' or name == 'flux' or name == 'spec' or name == 'sind'):
                            continue
                        if not name in liststrgfeatodim[l]:
                            nametotl = name + gdat.listnamerefr[q]
                            if name == 'etag':
                                continue
                            liststrgfeatodim[l].append(nametotl)
                            
                            if gdat.refrnumbelempopl[q] == 0:
                                continue

                            gdat.refrliststrgfeatonly[q][l].append(name)
                            if not nametotl in gdat.fittliststrgfeatextr[l]:
                                gdat.fittliststrgfeatextr[l].append(nametotl) 
                            #if name == 'reds':
                            #    for nametemp in ['lumi', 'dlos']:
                            #        nametemptemp = nametemp + gdat.listnamerefr[q]
                            #        if not nametemptemp in gdat.fittliststrgfeatextr[l]:
                            #            liststrgfeatodim[l].append(nametemp + gdat.listnamerefr[q])
                            #            gdat.fittliststrgfeatextr[l].append(nametemptemp)
        
        if gdat.exprtype == 'chan' and gdat.datatype == 'inpt':
            for l in indxpopl:
                if elemtype[l] == 'lghtpnts':
                    gdat.fittliststrgfeatextr[l].append('lumiwo08')
                    liststrgfeatodim[l].append('lumiwo08')

        # defaults
        liststrgpdfnmodu = [[] for l in indxpopl]
        liststrgfeatmodu = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght'): 
                if gdat.exprtype == 'ferm' and gdat.lgalcntr == 0.:
                    if l == 1:
                        liststrgpdfnmodu[l] += ['tmplnfwp']
                        liststrgfeatmodu[l] += ['lgalbgal']
                    if l == 2:
                        liststrgpdfnmodu[l] += ['tmplnfwp']
                        liststrgfeatmodu[l] += ['lgalbgal']
        
        # bring together two kinds of priors
        liststrgpdfnprio = [[] for l in indxpopl]
        liststrgfeatprio = [[] for l in indxpopl]
        for l in indxpopl:
            for strgfeat, strgpdfn in zip(liststrgfeatmodu[l], liststrgpdfnmodu[l]): 
                if not strgfeat in liststrgfeatprio[l]:
                    liststrgfeatprio[l].append(strgcomp)
                    liststrgpdfnprio[l].append(scalcomp)
            for strgcomp, scalcomp in zip(liststrgcomp[l], listscalcomp[l]):
                if not strgcomp in liststrgfeatprio[l]:
                    liststrgfeatprio[l].append(strgcomp)
                    liststrgpdfnprio[l].append(scalcomp)
        
        liststrgfeat = [[] for l in indxpopl]
        for l in indxpopl:
            for liststrg in [liststrgfeatprio[l], liststrgfeatodim[l], liststrgfeatsign[l]]:
                for strgthis in liststrg:
                    if not strgthis in liststrgfeat[l]:
                        liststrgfeat[l].append(strgthis)
        
        for l in indxpopl:
            for liststrgtemp in [liststrgfeat, liststrgfeatodim]:
                for strgfeat in liststrgtemp[l]:
                    if liststrgtemp[l].count(strgfeat) != 1:
                        print 'liststrgtemp'
                        print liststrgtemp
                        print
                        raise Exception('')

        # temp
        for l in indxpopl:
            if elemtype[l].startswith('lghtline'):
                liststrgfeat[l] += ['spec']
            if elemtype[l].startswith('lght'):
                liststrgfeat[l] += ['spec', 'specplot']
            if elemtype[l] == 'lens':
                liststrgfeat[l] += ['deflprof']
        
        liststrgfeateval = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('clus'):
                liststrgfeateval[l] = ['lgal', 'bgal', 'nobj']
            if elemtype[l] == 'clusvari':
                liststrgfeateval[l] += ['gwdt']
            if elemtype[l] == 'lens':
                liststrgfeateval[l] = ['lgal', 'bgal', 'defs', 'asca', 'acut']
            if elemtype[l].startswith('lghtline'):
                liststrgfeateval[l] = ['elin', 'spec']
            elif elemtype[l] == 'lghtgausbgrd':
                liststrgfeateval[l] = ['lgal', 'bgal', 'gwdt', 'spec']
            elif elemtype[l].startswith('lght'):
                liststrgfeateval[l] = ['lgal', 'bgal', 'spec']
    
        ## legends
        legdpopl = [[] for l in indxpopl]
        for l in indxpopl:
            if gdat.numbgrid > 1:
                if elemtype[l] == 'lghtpnts':
                    legdpopl[l] = 'FPS'
                if elemtype[l] == 'lghtgausbgrd':
                    legdpopl[l] = 'BGS'
            else:
                if elemtype[l] == 'lghtpntspuls':
                    legdpopl[l] = 'Pulsar'
                elif elemtype[l].startswith('lghtpntsagnn'):
                    legdpopl[l] = 'AGN'
                elif elemtype[l].startswith('lghtpnts'):
                    legdpopl[l] = 'PS'
            if elemtype[l] == 'lens':
                legdpopl[l] = 'Subhalo'
            if elemtype[l].startswith('clus'):
                legdpopl[l] = 'Cluster'
            if elemtype[l].startswith('lghtline'):
                legdpopl[l]= 'Line'
        
        for l in indxpopl:
            if 'gang' in liststrgcomp[l] and elemtype[l] != 'lghtpntspuls':
                setp_varbvalu(gdat, 'maxmgang', gdat.maxmgangdata, strgmodl=strgmodl)
        
        if strgmodl == 'true':
            setp_varbvalu(gdat, 'legdpopl', legdpopl, strgmodl='refr')
            legdpopl = gdat.refrlegdpopl
        else:
            setp_varbvalu(gdat, 'legdpopl', legdpopl, strgmodl=strgmodl)
            legdpopl = getattr(gdat, strgmodl + 'legdpopl')

        if strgmodl == 'true':
            indxpoplassc = [[] for l in indxpopl]
            for l in indxpopl:
                if numbpopl == 3 and elemtype[1] == 'lens':
                    indxpoplassc[l] = [l]
                else:
                    indxpoplassc[l] = indxpopl

        # variables for which two dimensional histograms will be calculated
        liststrgfeatcorr = [[] for l in indxpopl]
        if gdat.plotelemcorr:
            for l in indxpopl:
                for strgfeat in liststrgfeatodim[l]:
                    liststrgfeatcorr[l].append(strgfeat)
        
        # number of element parameters
        if numbpopl > 0:
            numbcomp = np.zeros(numbpopl, dtype=int)
            for l in indxpopl:
                numbcomp[l] = len(liststrgcomp[l])
            maxmnumbcomp = np.amax(numbcomp)
    
        # size of the auxiliary variable propobability density vector
        if maxmnumbelemtotl > 0:
            numblpri = 3 + maxmnumbcomp * numbpopl
        else:
            numblpri = 0
        if gdat.penalpridiff:
            numblpri += 1
        indxlpri = np.arange(numblpri)

        indxcomp = []
        for l in indxpopl:
            indxcomp.append(np.arange(numbcomp[l]))

        boolcompposi = [[] for l in indxpopl]
        for l in indxpopl:
            boolcompposi[l] = np.zeros(numbcomp[l], dtype=bool)
            if elemtype[l].startswith('lghtline'):
                boolcompposi[l][0] = True
            else:
                boolcompposi[l][0] = True
                boolcompposi[l][1] = True
        
        # number of transdimensional parameters
        numbtrappopl = np.zeros(numbpopl, dtype=int)
        numbtrappoplcuml = np.zeros(numbpopl, dtype=int)
        numbtrappoplcumr = np.zeros(numbpopl, dtype=int)
        for l in indxpopl:
            numbtrappopl[l] = maxmnumbelem[l] * numbcomp[l]
        for l in indxpopl:
            numbtrappoplcuml[l] = np.sum(numbtrappopl[:l])
            numbtrappoplcumr[l] = np.sum(numbtrappopl[:l+1])
        
        numbtrap = 0
        for l in indxpopl:
            numbtrap += np.sum(numbtrappopl[l])
        
        # list of strings across all populations
        ## element parameters
        liststrgcomptotl = retr_listconc(liststrgcomp)
        numbcomptotl = len(liststrgcomptotl)
        indxcomptotl = np.arange(numbcomptotl)
    
        ## element parameters
        liststrgfeatpriototl = retr_listconc(liststrgfeatprio)
        liststrgfeattotl = retr_listconc(liststrgfeat)
        liststrgfeatcorrtotl = retr_listconc(liststrgfeatcorr)
        
        ## one dimensional element features
        liststrgfeatodimtotl = retr_listconc(liststrgfeatodim)

        numbelemfeat = len(liststrgfeattotl)
        
        numbdeflsubhplot = 2
        numbdeflsingplot = numbdeflsubhplot
        if numbtrap > 0:
            numbdeflsingplot += 3

        listnamediff = []
        for c in indxback:
            listnamediff += ['back%04d' % c]
        if hostemistype != 'none':
            for e in indxsersfgrd:
                listnamediff += ['hostisf%d' % e]
        if lensmodltype != 'none':
            listnamediff += ['lens']
        
        listnameecom = deepcopy(listnamediff)
        for l in indxpopl:
            if boolelemsbrt[l]:
                if strgmodl == 'true' and numbelempopl[l] > 0 or strgmodl == 'fitt' and maxmnumbelempopl[l] > 0:
                    if not 'dfnc' in listnameecom:
                        listnameecom += ['dfnc']
                    if not 'dfncsubt' in listnameecom:
                        listnameecom += ['dfncsubt']
        listnameecomtotl = listnameecom + ['modl']
        
        listnamegcom = deepcopy(listnameecomtotl)
        if lensmodltype != 'none':
            listnamegcom += ['bgrd']
            if numbtrap > 0 and boolelemsbrtextsbgrdanyy:
                listnamegcom += ['bgrdgalx', 'bgrdexts']
        
        numbdiff = len(listnamediff)
        convdiff = np.zeros(numbdiff, dtype=bool)
        for k, namediff in enumerate(listnamediff):
            if not (gdat.boolthindata or psfnevaltype == 'none' or psfnevaltype == 'kern'):
                if namediff.startswith('back'):
                    indx = int(namediff[-4:])
                    convdiff[k] = not unifback[indx] 
                else:
                    convdiff[k] = True
        convdiffanyy = True in convdiff

        cntr = tdpy.util.cntr()
        
        liststrgfeatdefa = deepcopy(liststrgfeattotl)
        if boolelemlghtanyy:
            for strgfeat in ['sind', 'curv', 'expc'] + ['sindcolr%04d' % i for i in gdat.indxenerinde]:
                if not strgfeat in liststrgfeatdefa:
                    liststrgfeatdefa.append(strgfeat)

        if lensmodltype != 'none':
            redshost = getattr(gdat, strgmodl + 'redshost')
            redssour = getattr(gdat, strgmodl + 'redssour')
            adishost = gdat.adisobjt(redshost)
            adissour = gdat.adisobjt(redssour)
            adishostsour = adissour - (1. + redshost) / (1. + redssour) * adishost
            massfrombein = retr_massfrombein(gdat, adissour, adishost, adishostsour)
            mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
        
        if numbtrap > 0:

            # number of elements
            for l in indxpopl:
                dicttemp['indxfixpnumbelempop%d' % l] = cntr.incr()
            
            # hyperparameters
            ## mean number of elements
            for l in indxpopl:
                if maxmnumbelempopl[l] > 0:
                    dicttemp['indxfixpmeanelempop%d' % l] = cntr.incr()
            
            for strgfeat in liststrgfeatpriototl:
                try:
                    disttype = getattr(gdat, strgmodl + strgfeat + 'disttype')
                except:
                    if strgfeat == 'lgal' or strgfeat == 'bgal' or strgfeat == 'elin':
                        disttype = 'self'
                    try:
                        setattr(gdat, strgmodl + strgfeat + 'disttype', disttype)
                    except:
                        pass

            ## distribution shapes
            liststrgvarb = []
            for l in indxpopl:
                if maxmnumbelempopl[l] > 0:
                    for pdfnfeat, strgfeat in zip(liststrgpdfnprio[l], liststrgfeatprio[l]):
                        if pdfnfeat == 'exposcal' or pdfnfeat == 'dexpscal':
                            liststrgvarb += [strgfeat + 'distscal']
                        if pdfnfeat == 'powrslop':
                            liststrgvarb += [strgfeat + 'distslop']
                        if pdfnfeat == 'dpowslopbrek':
                            liststrgvarb += [strgfeat + 'distbrek']
                            liststrgvarb += [strgfeat + 'distsloplowr']
                            liststrgvarb += [strgfeat + 'distslopuppr']
                        if pdfnfeat == 'gausmean' or pdfnfeat == 'lnormean':
                            liststrgvarb += [strgfeat + 'distmean']
                        if pdfnfeat == 'gausstdv' or pdfnfeat == 'lnorstdv':
                            liststrgvarb += [strgfeat + 'diststdv']
                        if pdfnfeat == 'gausmeanstdv' or pdfnfeat == 'lnormeanstdv':
                            liststrgvarb += [strgfeat + 'distmean', strgfeat + 'diststdv']
            # temp
            for strgvarb in liststrgvarb:
                strgtemp = 'indxfixp' + strgvarb
                temp = np.zeros(numbpopl, dtype=int) - 1
                dicttemp[strgtemp] = temp

            for l in indxpopl:
                if maxmnumbelempopl[l] > 0:
                    for k, strgfeatprio in enumerate(liststrgfeatprio[l]):
                        if liststrgpdfnprio[l][k] == 'exposcal' or liststrgpdfnprio[l][k] == 'dexpscal':
                            dicttemp['indxfixp' + strgfeatprio + 'distscalpop%d' % l] = cntr.incr()
                            dicttemp['indxfixp' + strgfeatprio + 'distscal'][l] = dicttemp['indxfixp' + strgfeatprio + 'distscalpop%d' % l]
                        if liststrgpdfnprio[l][k] == 'gaumcons':
                            dicttemp['indxfixplgalbgaldistconspop%d' % l] = cntr.incr()
                            dicttemp['indxfixplgalbgaldistcons'][l] = dicttemp['indxfixplgalbgaldistconspop%d' % l]
                        if liststrgpdfnprio[l][k] == 'powrslop':
                            dicttemp['indxfixp' + strgfeatprio + 'distsloppop%d' % l] = cntr.incr()
                            dicttemp['indxfixp' + strgfeatprio + 'distslop'][l] = dicttemp['indxfixp' + strgfeatprio + 'distsloppop%d' % l]
                        if liststrgpdfnprio[l][k] == 'dpowslopbrek':
                            for nametemp in ['brek', 'sloplowr', 'slopuppr']:
                                dicttemp['indxfixp' + strgfeatprio + 'dist%spop%d' % (nametemp, l)] = cntr.incr()
                                dicttemp['indxfixp' + strgfeatprio + 'dist%s' % nametemp][l] = dicttemp['indxfixp' + strgfeatprio + 'dist%spop%d' % (nametemp, l)]
                        if liststrgpdfnprio[l][k] == 'gausmean' or liststrgpdfnprio[l][k] == 'gausmeanstdv' or \
                                        liststrgpdfnprio[l][k] == 'lnormean' or liststrgpdfnprio[l][k] == 'lnormeanstdv':
                            dicttemp['indxfixp' + strgfeatprio + 'distmeanpop%d' % l] = cntr.incr()
                            dicttemp['indxfixp' + strgfeatprio + 'distmean'][l] = dicttemp['indxfixp' + strgfeatprio + 'distmeanpop%d' % l]
                        if liststrgpdfnprio[l][k] == 'gausstdv' or liststrgpdfnprio[l][k] == 'gausmeanstdv' or \
                                        liststrgpdfnprio[l][k] == 'lnorstdv' or liststrgpdfnprio[l][k] == 'lnormeanstdv':
                            dicttemp['indxfixp' + strgfeatprio + 'diststdvpop%d' % l] = cntr.incr()
                            dicttemp['indxfixp' + strgfeatprio + 'diststdv'][l] = dicttemp['indxfixp' + strgfeatprio + 'diststdvpop%d' % l]
        
            dicttemp['indxfixpnumbelem'] = [[] for l in indxpopl]
            dicttemp['indxfixpmeanelem'] = []
            for strg, valu in dicttemp.iteritems():
                if strg[8:].startswith('numbelemp'):
                    indxpopltemp = int(strg[8:][11])
                    dicttemp['indxfixpnumbelem'][indxpopltemp] = valu
                if strg[8:].startswith('meanelemp'):
                    dicttemp['indxfixpmeanelem'].append(valu)
            dicttemp['indxfixpmeanelem'] = np.array(dicttemp['indxfixpmeanelem'])
            
            dicttemp['indxfixpdist'] = []
            for strg, valu in dicttemp.iteritems():
                if strg[12:16] == 'dist' and np.isscalar(valu):
                    dicttemp['indxfixpdist'].append(valu)
                
            dicttemp['indxfixpdist'] = np.array(dicttemp['indxfixpdist']) 
            dicttemp['indxfixphypr'] = np.array(list(dicttemp['indxfixpdist']) + list(dicttemp['indxfixpmeanelem']))
        
        if psfnevaltype != 'none':
            for m in gdat.indxevtt:
                for i in gdat.indxener:
                    dicttemp['indxfixpsigcen%02devt%d' % (i, m)] = cntr.incr()
                    if psfntype == 'doubking' or psfntype == 'singking':
                        dicttemp['indxfixpgamcen%02devt%d' % (i, m)] = cntr.incr()
                        if psfntype == 'doubking':
                            dicttemp['indxfixpsigten%02devt%d' % (i, m)] = cntr.incr()
                            dicttemp['indxfixpgamten%02devt%d' % (i, m)] = cntr.incr()
                            dicttemp['indxfixppsffen%02devt%d' % (i, m)] = cntr.incr()
        
            dicttemp['indxfixppsfp'] = []
            for strg, valu in dicttemp.iteritems():
                if strg.startswith('indxfixpsigce') \
                                            or strg.startswith('indxfixpsigte') or strg.startswith('indxfixpgamce') or strg.startswith('indxfixpgamte') or \
                                            strg.startswith('indxfixppsffe'):
                    dicttemp['indxfixppsfp'].append(valu)
            dicttemp['indxfixppsfp'] = np.array(dicttemp['indxfixppsfp']) 

            numbpsfptotlevtt = gdat.numbevtt * numbpsfptotl
            numbpsfptotlener = gdat.numbener * numbpsfptotl
            numbpsfp = numbpsfptotl * gdat.numbener * gdat.numbevtt
            indxpsfpform = np.arange(numbpsfpform)
            indxpsfptotl = np.arange(numbpsfptotl)
   
            indxpsfp = np.arange(numbpsfp)
            indxfixppsfp = np.sort(dicttemp['indxfixppsfp'])
            indxfixppsfpinit = indxfixppsfp[0]
        
        dicttemp['indxfixpbacp'] = []
        for c in indxback:
            if specback[c]:
                indx = cntr.incr()
                dicttemp['indxfixpbacpback%04d' % c] = indx
                dicttemp['indxfixpbacp'].append(indx)
            else:
                for i in gdat.indxener:
                    indx = cntr.incr()
                    dicttemp['indxfixpbacpback%04den%02d' % (c, i)] = indx
                    dicttemp['indxfixpbacp'].append(indx)
                    
        dicttemp['indxfixpbacp'] = np.array(dicttemp['indxfixpbacp'])
        
        # temp
        #dicttemp['indxfixpanglsour'] = []
        #dicttemp['indxfixpanglhost'] = []
        #dicttemp['indxfixpangllens'] = []
        
        if hostemistype != 'none':
            dicttemp['indxfixpspecsour'] = []
            dicttemp['indxfixpspechost'] = []

        if lensmodltype != 'none':
            dicttemp['indxfixplgalsour'] = cntr.incr()
            dicttemp['indxfixpbgalsour'] = cntr.incr()
            dicttemp['indxfixpfluxsour'] = cntr.incr()
            if gdat.numbener > 1:
                dicttemp['indxfixpsindsour'] = cntr.incr()
            dicttemp['indxfixpsizesour'] = cntr.incr()
            dicttemp['indxfixpellpsour'] = cntr.incr()
            dicttemp['indxfixpanglsour'] = cntr.incr()
        if hostemistype != 'none' or lensmodltype != 'none':
            for e in indxsersfgrd: 
                if hostemistype != 'none':
                    dicttemp['indxfixplgalhostisf%d' % e] = cntr.incr()
                    dicttemp['indxfixpbgalhostisf%d' % e] = cntr.incr()
                    dicttemp['indxfixpfluxhostisf%d' % e] = cntr.incr()
                    if gdat.numbener > 1:
                        dicttemp['indxfixpsindhostisf%d' % e] = cntr.incr()
                    dicttemp['indxfixpsizehostisf%d' % e] = cntr.incr()
                if lensmodltype != 'none':
                    dicttemp['indxfixpbeinhostisf%d' % e] = cntr.incr()
                if hostemistype != 'none':
                    dicttemp['indxfixpellphostisf%d' % e] = cntr.incr()
                    dicttemp['indxfixpanglhostisf%d' % e] = cntr.incr()
                    dicttemp['indxfixpserihostisf%d' % e] = cntr.incr()
        if lensmodltype != 'none':
            dicttemp['indxfixpsherextr'] = cntr.incr()
            dicttemp['indxfixpsangextr'] = cntr.incr()
            dicttemp['indxfixpsour'] = []
        
        if lensmodltype != 'none' and hostemistype == 'none':
            raise Exception('Lensing cannot be modeled without host galaxy emission.')

        liststrgcomplens = ['hostlght', 'hostlens', 'sour', 'extr']
        dictname = dict()
        for strgcomplens in liststrgcomplens:
            dictname['liststrg' + strgcomplens] = []
            dicttemp['indxfixp' + strgcomplens] = []
        if lensmodltype != 'none' or hostemistype != 'none':
            dictname['liststrghostlght'] += ['lgalhost', 'bgalhost', 'ellphost', 'anglhost']
            dictname['liststrghostlens'] += ['lgalhost', 'bgalhost', 'ellphost', 'anglhost']
        if hostemistype != 'none':
            dictname['liststrghostlght'] += ['fluxhost', 'sizehost', 'serihost']
            if gdat.numbener > 1:
                dictname['liststrghostlght'] += ['sindhost']
        if lensmodltype != 'none':
            dictname['liststrghostlens'] += ['beinhost']
            dictname['liststrgextr'] += ['sherextr', 'sangextr']
            dictname['liststrgsour'] += ['lgalsour', 'bgalsour', 'fluxsour', 'sizesour', 'ellpsour', 'anglsour']
            if gdat.numbener > 1:
                dictname['liststrgsour'] += ['sindsour']

        for strg, valu in dicttemp.iteritems():
            
            if isinstance(valu, list) or isinstance(valu, np.ndarray):
                continue
            
            for strgcomplens in liststrgcomplens:
                setp_indxfixparry(dicttemp, strgcomplens, strg, dictname['liststrg' + strgcomplens], valu)
                
            if strg[8:].startswith('fluxsour') or strg[8:].startswith('sindsour'):
                dicttemp['indxfixpspecsour'].append(valu)

            if strg[8:].startswith('fluxhost') or strg[8:].startswith('sindhost'):
                dicttemp['indxfixpspechost'].append(valu)
            
        dicttemp['indxfixphost'] = dicttemp['indxfixphostlght'] + dicttemp['indxfixphostlens']
        dicttemp['indxfixplenp'] = dicttemp['indxfixphost'] + dicttemp['indxfixpsour'] + dicttemp['indxfixpextr']

        # number of fixed-dimension parameters
        numbfixp = cntr.incr(0)
        
        # indices of fixed-dimension parameters
        indxfixp = np.arange(numbfixp)
       
        indxfixpstdv = indxfixp[numbpopl:]

        ## number of model spectral parameters for each population
        #numbspep = np.empty(numbpopl, dtype=int)
        #liststrgspep = [[] for l in range(numbpopl)]
        #for l in indxpopl:
        #    if gdat.numbener > 1:
        #        liststrgspep[l] += ['sind']
        #        if spectype[l] == 'expc':
        #            liststrgspep[l] += ['expc']
        #        if spectype[l] == 'curv':
        #            liststrgspep[l] = ['curv']
        #    numbspep[l] = len(liststrgspep[l]) 

        # total number of parameters
        numbpara = numbfixp + numbtrap
        
        indxsamptrap = np.arange(numbfixp, numbpara)
        if numbtrap > 0:
            indxsampcompinit = indxsamptrap[0]
        else:
            indxsampcompinit = -1
        indxpara = np.arange(numbpara)

    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr' and len(strg) % 4 == 0:
            if hasattr(gdat, strgmodl + strg):
                pass
                #print 'Already has %s' % (strgmodl + strg)
                #print
                #print
            else:
                setattr(gdat, strgmodl + strg, valu)

    for strg, valu in dicttemp.iteritems():
        # temp
        if isinstance(valu, np.ndarray) and strg[12:16] != 'dist':
            valu = np.sort(valu)
        setattr(gdat, strgmodl + strg, valu)


def setp_indxfixparry(dicttemp, name, strg, liststrg, valu):
    
    for strgtemp in liststrg:
        if strg[8:].startswith(strgtemp):
            if isinstance(valu, list):
                for valutemp in valu:
                    dicttemp['indxfixp' + name].append(valutemp)
            else:
                dicttemp['indxfixp' + name].append(valu)
           

def plot_lens(gdat):
    
    if gdat.fittboolelemdeflsubh:
        xdat = gdat.binsangl[1:] * gdat.anglfact
        lablxdat = gdat.lablgangtotl
        
        listdeflscal = np.array([4e-2, 4e-2, 4e-2]) / gdat.anglfact
        listanglscal = np.array([0.05, 0.1, 0.05]) / gdat.anglfact
        listanglcutf = np.array([1.,    1.,  10.]) / gdat.anglfact
        listasym = [False, False, False]
        listydat = []
        for deflscal, anglscal, anglcutf, asym in zip(listdeflscal, listanglscal, listanglcutf, listasym):
            listydat.append(retr_deflcutf(gdat.binsangl[1:], deflscal, anglscal, anglcutf, asym=asym) * gdat.anglfact)
        
        for scalxdat in ['self', 'logt']:
            path = gdat.pathinitintr + 'deflcutf' + scalxdat + '.pdf'
            tdpy.util.plot_gene(path, xdat, listydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, \
                                                             lablydat=r'$\alpha_n$ [$^{\prime\prime}$]', limtydat=[1e-3, 1.5e-2], limtxdat=[None, 2.])
       
        # pixel-convoltuion of the Sersic profile
        # temp -- y axis labels are wrong, should be per solid angle
        xdat = gdat.binslgalsers * gdat.anglfact
        for n in range(gdat.numbindxsers + 1):
            for k in range(gdat.numbhalfsers + 1):
                if k != 5:
                    continue
                path = gdat.pathinitintr + 'sersprofconv%04d%04d.pdf' % (n, k)
                tdpy.util.plot_gene(path, xdat, gdat.sersprof[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                #path = gdat.pathinitintr + 'sersprofcntr%04d%04d.pdf' % (n, k)
                #tdpy.util.plot_gene(path, xdat, gdat.sersprofcntr[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.pdf' % (n, k)
                tdpy.util.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.pdf' % (n, k)
                tdpy.util.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], scalxdat='logt', \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
       
        xdat = gdat.binsangl * gdat.anglfact
        listspec = np.array([1e-19, 1e-18, 1e-18, 1e-18]) / gdat.anglfact
        listsize = np.array([0.3, 1., 1., 1.]) / gdat.anglfact
        listindx = np.array([4., 2., 4., 10.])
        listydat = []
        listlegd = []
        for spec, size, indx in zip(listspec, listsize, listindx):
            listydat.append(spec * retr_sbrtsersnorm(gdat.binsangl, size, indxsers=indx))
            listlegd.append('$R_e = %.3g ^{\prime\prime}, n = %.2g$' % (size * gdat.anglfact, indx))
        path = gdat.pathinitintr + 'sersprof.pdf'
        tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, \
                                                                                                   listlegd=listlegd, listhlin=1e-7, limtydat=[1e-8, 1e0])
    
        minmredshost = 0.01
        maxmredshost = 0.4
        minmredssour = 0.01
        maxmredssour = 2.
        numbreds = 200
        retr_axis(gdat, 'redshost', minmredshost, maxmredshost, numbreds)
        retr_axis(gdat, 'redssour', minmredssour, maxmredssour, numbreds)
        
        gdat.meanadishost = np.empty(numbreds)
        for k in range(numbreds):
            gdat.meanadishost[k] = gdat.adisobjt(gdat.meanredshost[k])
        
        asca = 0.1 / gdat.anglfact
        acut = 1. / gdat.anglfact
    
        minmmass = np.zeros((numbreds + 1, numbreds + 1))
        maxmmass = np.zeros((numbreds + 1, numbreds + 1))
        for k, redshost in enumerate(gdat.binsredshost):
            for n, redssour in enumerate(gdat.binsredssour):
                if redssour > redshost:
                    adishost = gdat.adisobjt(redshost)
                    adissour = gdat.adisobjt(redssour)
                    adishostsour = adissour - (1. + redshost) / (1. + redssour) * adishost
                    factmcutfromdefs = retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut)
                    minmmass[n, k] = np.log10(factmcutfromdefs * gdat.minmdefs)
                    maxmmass[n, k] = np.log10(factmcutfromdefs * gdat.maxmdefs)
       
        #valulevl = np.linspace(7.5, 9., 5)
        valulevl = [7.0, 7.3, 7.7, 8., 8.6]
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        cont = axis.contour(gdat.binsredshost, gdat.binsredssour, minmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=20, fmt='%.3g')
        axis.set_xlabel(r'$z_{\rm{hst}}$')
        axis.set_ylabel(r'$z_{\rm{src}}$')
        axis.set_title(r'$M_{c,min}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsminm.pdf'
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
        
        valulevl = np.linspace(9., 11., 20)
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.imshow(maxmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=9., vmax=11.)
        cont = axis.contour(gdat.binsredshost, gdat.binsredssour, maxmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=15, fmt='%.3g')
        axis.set_xlabel('$z_{hst}$')
        axis.set_ylabel('$z_{src}$')
        axis.set_title(r'$M_{c,max}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsmaxm.pdf'
        plt.colorbar(imag) 
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.plot(gdat.meanredshost, gdat.meanadishost * gdat.sizepixl * 1e-3)
        axis.plot(gdat.meanredshost, gdat.meanadishost * 2. * gdat.maxmgangdata * 1e-3)
        axis.set_xlabel('$z_h$')
        axis.set_yscale('log')
        axis.set_ylabel(r'$\lambda$ [kpc]')
        path = gdat.pathinitintr + 'wlenreds.pdf'
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        fracacutasca = np.logspace(-1., 2., 20)
        mcut = retr_mcutfrommscl(fracacutasca)
        axis.lognp.log(fracacutasca, mcut)
        axis.set_xlabel(r'$\tau_n$')
        axis.set_ylabel(r'$M_{c,n} / M_{0,n}$')
        axis.axhline(1., ls='--')
        path = gdat.pathinitintr + 'mcut.pdf'
        plt.tight_layout()
        figr.savefig(path)
        plt.close(figr)
       

def retr_listrtagprev(strgcnfg):
    
    # list of PCAT run plot outputs
    pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
    listrtag = fnmatch.filter(os.listdir(pathimag), '2*')
    
    listrtagprev = []
    for rtag in listrtag:
        strgstat = os.environ["PCAT_DATA_PATH"] + '/data/outp/' + rtag
        
        if chec_statfile(rtag, 'gdatmodi', 'post', verbtype=0) and strgcnfg + '_' + rtag[16:].split('_')[-1] == rtag[16:]:
            listrtagprev.append(rtag) 
    
    listrtagprev.sort()

    return listrtagprev


def setp_fixp(gdat, strgmodl='fitt'):
    
    if gdat.verbtype > 0:
        if strgmodl == 'true':
            strgtemp = 'true'
        if strgmodl == 'fitt':
            strgtemp = 'fitting'
        print 'Building the %s model...' % strgtemp
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    
    psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')
    
    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbpopl > 0:
        elemtype = getattr(gdat, strgmodl + 'elemtype')
        numbcomp = getattr(gdat, strgmodl + 'numbcomp')
        numbtrappoplcuml = getattr(gdat, strgmodl + 'numbtrappoplcuml')
        legdpopl = getattr(gdat, strgmodl + 'legdpopl')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    numbfixp = getattr(gdat, strgmodl + 'numbfixp')
    numbback = getattr(gdat, strgmodl + 'numbback')
    numbpara = getattr(gdat, strgmodl + 'numbpara')
    listnameback = getattr(gdat, strgmodl + 'listnameback')
    boolelemlght = getattr(gdat, strgmodl + 'boolelemlght')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    boolelemsbrt = getattr(gdat, strgmodl + 'boolelemsbrt')
    boolelemsbrtpnts = getattr(gdat, strgmodl + 'boolelemsbrtpnts')
    boolelemsbrtdfnc = getattr(gdat, strgmodl + 'boolelemsbrtdfnc')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    lablelemextn = getattr(gdat, strgmodl + 'lablelemextn')
    if numbtrap > 0:
        namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
    if lensmodltype != 'none' or hostemistype != 'none':
        indxsersfgrd = getattr(gdat, strgmodl + 'indxsersfgrd')
    
    listlegdback = []
    listlablback = []
    for nameback in listnameback:
        if nameback == 'isot':
            listlegdback.append('Isotropic')
            listlablback.append(r'$\mathcal{I}$')
    	if nameback == 'fdfm':
    	    listlegdback.append('FDM')
            listlablback.append(r'$\mathcal{D}$')
    	if nameback == 'dark':
    	    listlegdback.append('NFW')
            listlablback.append(r'$\mathcal{D}_{dark}$')
    	if nameback == 'part':
    	    listlegdback.append('Particle Back.')
            listlablback.append(r'$\mathcal{I}_p$')
    
    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    psfntype = getattr(gdat, strgmodl + 'psfntype')
    if psfnevaltype != 'none':
        numbpsfptotl = getattr(gdat, strgmodl + 'numbpsfptotl')
    specback = getattr(gdat, strgmodl + 'specback')
    indxbackbacp = getattr(gdat, strgmodl + 'indxbackbacp')
    indxenerbacp = getattr(gdat, strgmodl + 'indxenerbacp')

    liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    scalmeanelem = getattr(gdat, strgmodl + 'scalmeanelem')

    # construct the fixed parameter structure
    namefixp = np.zeros(numbfixp, dtype=object)
    lablfixp = np.zeros(numbfixp, dtype=object)
    lablfixpunit = np.zeros(numbfixp, dtype=object)
    lablfixpunit[:] = ''
    scalfixp = np.zeros(numbfixp, dtype=object)
    minmfixp = np.zeros(numbfixp)
    maxmfixp = np.zeros(numbfixp)
    factfixp = np.zeros(numbfixp)
    meanfixp = np.zeros(numbfixp)
    stdvfixp = np.zeros(numbfixp)
    cdfnminmfixp = np.empty(numbfixp)
    cdfndifffixp = np.empty(numbfixp)
    factfixpplot = np.ones(numbfixp)
    
    for strg, k in gdat.__dict__.iteritems():
        
        strgvarb = strg[4:]

        if not isinstance(k, int) or not strgvarb.startswith('indxfixp') or strgvarb == 'indxfixppsfpinit' or strg[:4] != strgmodl:
            continue
        
        strgvarb = strgvarb[8:]
        namefixp[k] = strgvarb
        if 'pop' in strg:
            
            if numbpopl == 1:
                strgpopl = ''
                strgpoplcomm = ''
            else:
                strgpopl = '%s' % strg[-1]
                strgpoplcomm = ',%s' % strg[-1]
            indxpopltemp = int(strg.split('pop')[1][0])
            if namefixp[k].startswith('numbelem'):
                lablfixp[k] = '$N_{%s%s}$' % (lablelemextn[indxpopltemp], strgpoplcomm)
                scalfixp[k] = 'pois'
            if namefixp[k].startswith('meanelem'):
                lablfixp[k] = r'$\mu_{%s%s}$' % (lablelemextn[indxpopltemp], strgpoplcomm)
                scalfixp[k] = scalmeanelem
    
            if namefixp[k].startswith('gangdistslopp'):
                lablfixp[k] = r'$\beta_{\theta%s}$' % strgpoplcomm
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('spatdistconsp'):
                lablfixp[k] = r'$\gamma_{c%s}$' % strgpoplcomm
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('bgaldistscalp'):
                lablfixp[k] = r'$\gamma_{b%s}$' % strgpoplcomm
                scalfixp[k] = 'self'
                lablfixpunit[k] = gdat.lablgangunit
            
            if namefixp[k].startswith('gwdtdistslopp'):
                lablfixp[k] = r'$\beta_{%s%s}$' % (strgpoplcomm, gdat.lablgwdt)
                scalfixp[k] = 'logt'

            if namefixp[k][4:].startswith('distslopp') and namefixp[k][:4] == namefeatampl[indxpopltemp]:
                lablfixp[k] = r'$\beta_{%s}$' % strgpopl
                scalfixp[k] = getattr(gdat, strgmodl + 'scal' + namefixp[k])
    
            # normal and log-normal ditributed
            for nametemp in ['sind', 'curv', 'per0', 'magf', 'lumi']:
                if namefixp[k].startswith(nametemp + 'distmean'):
                    lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, getattr(gdat, 'labl' + nametemp))
                    scalfixp[k] = 'logt'
                
                if namefixp[k].startswith(nametemp + 'diststdv'):
                    lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, getattr(gdat, 'labl' + nametemp))
                    scalfixp[k] = 'logt'
            
            # power-law distributed
            for nametemp in ['dlos', 'dglc']:
                if namefixp[k].startswith(nametemp + 'distslop'):
                    lablfixp[k] = r'$\beta_{%s%s}$' % (strgpoplcomm, getattr(gdat, 'labl' + nametemp))
                    scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('fluxdistbrek'):
                lablfixp[k] = r'$f_{b,%s}$' % (strgpoplcomm)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('fluxdistsloplowr'):
                lablfixp[k] = r'$\beta_{f0,%s}$' % (strgpoplcomm)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('fluxdistslopuppr'):
                lablfixp[k] = r'$\beta_{f1,%s}$' % (strgpoplcomm)
                scalfixp[k] = 'logt'
            
            # double power-law distributed
            if namefixp[k].startswith('lum0distbrek'):
                lablfixp[k] = r'$L_{b,%s}$' % (strgpoplcomm)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('lum0distsloplowr'):
                lablfixp[k] = r'$\beta_{L0,%s%s}$' % (strgpoplcomm, gdat.labllum0)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('lum0distslopuppr'):
                lablfixp[k] = r'$\beta_{L1,%s%s}$' % (strgpoplcomm, gdat.labllum0)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('expcdistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablexpc)
                scalfixp[k] = 'logt'
                lablfixpunit[k] = gdat.lablenerunit
            
            if namefixp[k].startswith('expcdiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablexpc)
                scalfixp[k] = 'logt'
                lablfixpunit[k] = gdat.lablenerunit
            
            if namefixp[k].startswith('ascadistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablasca)
                scalfixp[k] = 'logt'
                factfixpplot[k] = gdat.anglfact
                lablfixpunit[k] = gdat.lablgangunit
            
            if namefixp[k].startswith('ascadiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablasca)
                scalfixp[k] = 'logt'
                factfixpplot[k] = gdat.anglfact
                lablfixpunit[k] = gdat.lablgangunit
            
            if namefixp[k].startswith('acutdistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablacut)
                scalfixp[k] = 'logt'
                factfixpplot[k] = gdat.anglfact
                lablfixpunit[k] = gdat.lablgangunit
            
            if namefixp[k].startswith('acutdiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablacut)
                scalfixp[k] = 'logt'
                factfixpplot[k] = gdat.anglfact
                lablfixpunit[k] = gdat.lablgangunit
            
        if strg[:-1].endswith('evt'):
            
            if gdat.psfninfoprio:
                scalfixp[k] = 'gaus'
                n = k - getattr(gdat, strgmodl + 'indxfixpsigcen00evt0')
                meanfixp[k] = gdat.psfpexpr[n]
                stdvfixp[k] = 0.1 * gdat.psfpexpr[n]
            else:
                if strgvarb.startswith('sig'):
                    scalfixp[k] = 'logt'
                if strgvarb.startswith('gam'):
                    scalfixp[k] = 'logt'
                if strgvarb.startswith('psff'):
                    scalfixp[k] = 'self'
            
            # strings for PSF parameters
            if strgvarb.startswith('sig'):
                strgvarbtemp = '\sigma'
                factfixpplot[k] = gdat.anglfact
                lablfixpunit[k] = gdat.lablgangunit
            if strgvarb.startswith('gam'):
                strgvarbtemp = '\gamma'
            if strgvarb.startswith('psff'):
                strgvarbtemp = 'f'
            if strgvarb.startswith('sig') and psfntype == 'doubgaus' or psfntype == 'gausking' or psfntype == 'doubking':
                strgcomptemp = 'c'
            elif strgvarb.startswith('gam') and psfntype == 'gausking' or psfntype == 'doubking':
                strgcomptemp = 't'
            else:
                strgcomptemp = ''
            if gdat.numbener > 1:
                indxenertemp = gdat.indxenerincl[((k - getattr(gdat, strgmodl + 'indxfixpsigcen00evt0')) % (gdat.numbener * numbpsfptotl)) // numbpsfptotl]
                strgenertemp = '%s' % indxenertemp
            else:
                strgenertemp = ''
            if gdat.numbevtt > 1:
                indxevtttemp = gdat.indxevttincl[(k - getattr(gdat, strgmodl + 'indxfixpsigcen00evt0')) // (gdat.numbener * numbpsfptotl)]
                strgevtttemp = '%s' % indxevtttemp
            else:
                strgevtttemp = ''
            lablfixp[k] = r'$%s^{%s}_{%s%s}$' % (strgvarbtemp, strgcomptemp, strgenertemp, strgevtttemp)
        
        if strgvarb.startswith('bacp'):
            
            # indices of the background
            indxbacktemp = int(strgvarb[8:12])
            
            name = 'bacpback%04d' % indxbacktemp
            
            if specback[indxbacktemp]:
                strgenertemp = ''
            else:
                indxenertemp = int(strgvarb[14:16])
                if gdat.numbener > 1:
                    strgenertemp = '%d' % indxenertemp
                else:
                    strgenertemp = ''
                name += 'en%02d' % indxenertemp
            
            if numbback > 1:
                strgbacktemp = '%d' % indxbacktemp
            else:
                strgbacktemp = ''
            
            lablfixp[k] = '$A_{%s%s}$' % (strgbacktemp, strgenertemp)
            lablfixpunit[k] = gdat.lablsbrtunit
            
            try:
                scalfixp[k] = getattr(gdat, strgmodl + 'scal' + name)
                if gdat.verbtype > 0:
                    print 'Received custom scaling for %s: %s' % (name, scalfixp[k])
            except:
                scalfixp[k] = 'logt'
            
        if lensmodltype == 'full' or lensmodltype == 'host':
            if k in getattr(gdat, strgmodl + 'indxfixplenp'):
                if strgvarb.startswith('lgalsour'):
                    lablfixp[k] = r'$\theta_{\rm{1,src}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('bgalsour'):
                    lablfixp[k] = r'$\theta_{\rm{2,src}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('fluxsour'):
                    lablfixp[k] = r'$f_{\rm{src}}$'
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if gdat.numbener > 1:
                    if strgvarb.startswith('sindsour'):
                        lablfixp[k] = r'$s_{\rm{src}}$'
                        scalfixp[k] = 'self'
                if strgvarb.startswith('sizesour'):
                    lablfixp[k] = r'$\theta_{\rm{e,src}}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('ellpsour'):
                    lablfixp[k] = r'$\epsilon_{\rm{src}}$'
                    scalfixp[k] = 'self'
                if strgvarb.startswith('anglsour'):
                    lablfixp[k] = r'$\phi_{\rm{src}}$'
                    scalfixp[k] = 'self'
                if strgvarb.startswith('lgalhost'):
                    lablfixp[k] = r'$\theta_{1,\rm{hst}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('bgalhost'):
                    lablfixp[k] = r'$\theta_{2,\rm{hst}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('fluxhost'):
                    lablfixp[k] = r'$f_{\rm{hst}}$'
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if gdat.numbener > 1:
                    if strgvarb.startswith('sindhost'):
                        lablfixp[k] = r'$s_{\rm{hst}}$'
                        scalfixp[k] = 'self'
                if strgvarb.startswith('sizehost'):
                    lablfixp[k] = r'$\theta_{\rm{e,hst}}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('beinhost'):
                    lablfixp[k] = r'$\theta_{\rm{E,hst}}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('ellphost'):
                    lablfixp[k] = r'$\epsilon_{\rm{hst}}$'
                    scalfixp[k] = 'self'
                if strgvarb.startswith('anglhost'):
                    lablfixp[k] = r'$\phi_{\rm{hst}}$'
                    scalfixp[k] = 'self'
                if strgvarb.startswith('serihost'):
                    lablfixp[k] = r'$n_{\rm{S,hst}}$'
                    scalfixp[k] = 'self'
                if strgvarb.startswith('sherextr'):
                    lablfixp[k] = r'$\gamma_{ext}$'
                    scalfixp[k] = 'self'
                if strgvarb.startswith('sangextr'):
                    lablfixp[k] = r'$\phi_{ext}$'
                    scalfixp[k] = 'self'
        
        if scalfixp[k] == 'pois' or scalfixp[k] == 'self' or scalfixp[k] == 'logt' or scalfixp[k] == 'atan':
            minmfixp[k] = getattr(gdat, strgmodl + 'minm' + strgvarb)
            maxmfixp[k] = getattr(gdat, strgmodl + 'maxm' + strgvarb)
        if scalfixp[k] == 'gaus' or scalfixp[k] == 'eerr':
            meanfixp[k] = getattr(gdat, strgmodl + 'mean' + strgvarb)
            stdvfixp[k] = getattr(gdat, strgmodl + 'stdv' + strgvarb)
        if scalfixp[k] == 'self':
            factfixp[k] = maxmfixp[k] - minmfixp[k]
        if scalfixp[k] == 'logt':
            factfixp[k] = np.log(maxmfixp[k] / minmfixp[k])
        if scalfixp[k] == 'atan':
            factfixp[k] = np.arctan(maxmfixp[k]) - np.arctan(minmfixp[k])
        if scalfixp[k] == 'gaus':
            minmfixp[k] = meanfixp[k] - 3. * stdvfixp[k]
            maxmfixp[k] = meanfixp[k] + 3. * stdvfixp[k]
        if scalfixp[k] == 'eerr':
            cdfnminmfixp[k], cdfndifffixp[k] = retr_eerrnorm(minmfixp[k], maxmfixp[k], meanfixp[k], stdvfixp[k])
        
        if not np.isfinite(factfixp[k]):
            print 'scalfixp[k]'
            print scalfixp[k]
            print 'namefixp[k]'
            print namefixp[k]
            print 'minmfixp[k]'
            print minmfixp[k]
            raise Exception('')
    
    listindxfixpscal = {}
    for scaltype in gdat.listscaltype:
        listindxfixpscal[scaltype] = np.where(scaltype == scalfixp)[0]

    # background templates
    listlegdsbrt = deepcopy(listlegdback)
    numblablsbrt = 0
    for l in indxpopl:
        if boolelemsbrt[l]:
            listlegdsbrt.append(legdpopl[l])
            listlegdsbrt.append(legdpopl[l] + ' subt')
            numblablsbrt += 2
    if lensmodltype != 'none':
        listlegdsbrt.append('Source')
        numblablsbrt += 1
    if hostemistype != 'none':
        for e in indxsersfgrd:
            listlegdsbrt.append('Host %d' % e)
            numblablsbrt += 1
    if numbpopl > 0:
        if 'clus' in elemtype or 'clusvari' in elemtype:
            listlegdsbrt.append('Uniform')
            numblablsbrt += 1
    
    listlegdsbrtspec = ['Data']
    listlegdsbrtspec += deepcopy(listlegdsbrt)
    if len(listlegdsbrt) > 1:
        listlegdsbrtspec.append('Total Model')
    
    numblablsbrtspec = len(listlegdsbrtspec)
    
    namepara = np.zeros(numbpara, dtype=object)
    lablpara = np.zeros(numbpara, dtype=object)
    scalpara = np.zeros(numbpara, dtype=object)
    factplotpara = np.zeros(numbpara)
    
    listfactplotcomp = []
    listlablcomp = []
    for l in indxpopl:
        listfactplotcomp.append(np.zeros(numbcomp[l]))
        listlablcomp.append(np.zeros(numbcomp[l], dtype=object))
        for k, namecomp in enumerate(liststrgcomp[l]):
            listfactplotcomp[l][k] = getattr(gdat, 'fact' + namecomp + 'plot')
            listlablcomp[l][k] = getattr(gdat, 'labl' + namecomp + 'totl')

    if strgmodl == 'fitt':
        for k in indxfixp:
            setattr(gdat, 'labl' + namefixp[k], lablfixp[k])
            setattr(gdat, 'fact' + namefixp[k] + 'plot', factfixpplot[k])
   
    namepara[indxfixp] = namefixp
    lablpara[indxfixp] = lablfixp
    scalpara[indxfixp] = scalfixp
    factplotpara[indxfixp] = factfixpplot
    for k in range(numbtrap):
        for l in indxpopl:  
            if k >= numbtrappoplcuml[l]:
                indxpopltemp = l
                indxelemtemp = (k - numbtrappoplcuml[indxpopltemp]) // numbcomp[indxpopltemp]
                indxcomptemp = (k - numbtrappoplcuml[indxpopltemp]) % numbcomp[indxpopltemp]
                break
        namepara[numbfixp+k] = '%spop%d%04d' % (liststrgcomp[indxpopltemp][indxcomptemp], indxpopltemp, indxelemtemp)
        lablpara[numbfixp+k] = listlablcomp[indxpopltemp][indxcomptemp]
        scalpara[numbfixp+k] = listscalcomp[indxpopltemp][indxcomptemp]
        factplotpara[numbfixp+k] = listfactplotcomp[indxpopltemp][indxcomptemp]
    
    #listscalparatotl = list(np.unique(scalpara))
    #listnameparatotl = list(np.unique(namepara))
    # limits for the Gaussian distributed parameters
    numbstdvgaus = 4.

    # determine minima and maxima for normal or log-normal distributed parameters
    for k, name in enumerate(namepara):
        if scalpara[k] == 'gaus' or scalpara[k] == 'lnormeanstdv':
            if k >= numbfixp:
                nametemp = name[:-12]
                for l in indxpopl:
                    if nametemp in liststrgcomp[l]:
                        minm = 1e100
                        maxm = -1e100
                        indx = liststrgcomp[l].index(nametemp)
                        if listscalcomp[l][indx] == 'gaus':
                            minmtemp = getattr(gdat, strgmodl + 'minm' + nametemp + 'distmeanpop%d' % l) - \
                                                    numbstdvgaus * getattr(gdat, strgmodl + 'maxm' + nametemp + 'diststdvpop%d' % l)
                            maxmtemp = getattr(gdat, strgmodl + 'maxm' + nametemp + 'distmeanpop%d' % l) + \
                                                    numbstdvgaus * getattr(gdat, strgmodl + 'maxm' + nametemp + 'diststdvpop%d' % l)
                            minm = min(minmtemp, minm)
                            maxm = max(maxmtemp, maxm)
                        if listscalcomp[l][indx] == 'lnormeanstdv':
                            minmtemp = np.exp(np.log(getattr(gdat, strgmodl + 'minm' + nametemp + 'distmeanpop%d' % l)) - \
                                                                        numbstdvgaus * getattr(gdat, strgmodl + 'maxm' + nametemp + 'diststdvpop%d' % l))
                            maxmtemp = np.exp(np.log(getattr(gdat, strgmodl + 'maxm' + nametemp + 'distmeanpop%d' % l)) + \
                                                                        numbstdvgaus * getattr(gdat, strgmodl + 'maxm' + nametemp + 'diststdvpop%d' % l))
                            minm = min(minmtemp, minm)
                            maxm = max(maxmtemp, maxm)

            else:
                nametemp = name
                minm = getattr(gdat, strgmodl + 'mean' + name) - numbstdvgaus * getattr(gdat, strgmodl + 'stdv' + name)
                maxm = getattr(gdat, strgmodl + 'mean' + name) + numbstdvgaus * getattr(gdat, strgmodl + 'stdv' + name)
            setattr(gdat, strgmodl + 'minm' + nametemp, minm)
            setattr(gdat, strgmodl + 'maxm' + nametemp, maxm)

    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr':
            if len(strg) % 4 == 0:
                setattr(gdat, strgmodl + strg, valu)


def make_legd(axis, offs=None, loca=1, numbcols=1, ptch=None, line=None):
   
    hand, labl = axis.get_legend_handles_labels()
    legd = axis.legend(hand, labl, fancybox=True, frameon=True, bbox_to_anchor=offs, bbox_transform=axis.transAxes, ncol=numbcols, loc=loca, labelspacing=1, handlelength=2)
    try:
        legd.get_frame().set_fill(True)
        legd.get_frame().set_facecolor('white')
    except:
        print 'Legend failed...'


def setp_namevarbsing(gdat, strgmodl, strgvarb, popl, ener, evtt, back, isfr):
    
    if popl == 'full':
        indxpopltemp = getattr(gdat, strgmodl + 'indxpopl')
    elif popl != 'none':
        indxpopltemp = [popl]
    
    if ener == 'full':
        indxenertemp = gdat.indxener
    elif ener != 'none':
        indxenertemp = [ener]
    
    if evtt == 'full':
        indxevtttemp = gdat.indxevtt
    elif evtt != 'none':
        indxevtttemp = [evtt]
    
    if back == 'full':
        indxbacktemp = getattr(gdat, strgmodl + 'indxback')
    elif isinstance(back, int):
        indxbacktemp = np.array([back])
    
    if isfr == 'full':
        indxisfrtemp = getattr(gdat, strgmodl + 'indxsersfgrd')
    
    liststrgvarb = []
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none' and isfr != 'none':
        for e in indxisfrtemp:
            liststrgvarb.append(strgvarb + 'isf%d' % e)
    if popl != 'none' and ener == 'none' and evtt == 'none' and back == 'none':
        for l in indxpopltemp:
            liststrgvarb.append(strgvarb + 'pop%d' % l)
    if popl != 'none' and ener == 'none' and evtt == 'none' and back == 'none':
        for l in indxpopltemp:
            liststrgvarb.append(strgvarb + 'pop%d' % l)
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none':
        liststrgvarb.append(strgvarb + '')
    if popl == 'none' and ener != 'none' and evtt != 'none' and back == 'none':
        for i in indxenertemp:
            for m in indxevtttemp:
                liststrgvarb.append(strgvarb + 'en%02devt%d' % (i, m))
    if popl == 'none' and ener != 'none' and evtt == 'none' and back != 'none':
        for c in indxbacktemp:
            for i in indxenertemp:
                liststrgvarb.append(strgvarb + 'back%04den%02d' % (c, i))
    if popl == 'none' and ener == 'none' and evtt == 'none' and back != 'none':
        for c in indxbacktemp:
            liststrgvarb.append(strgvarb + 'back%04d' % c)
    if popl == 'none' and ener != 'none' and evtt == 'none' and back != 'none':
        for c in indxbacktemp:
            for i in indxenertemp:
                liststrgvarb.append(strgvarb + 'back%04den%02d' % (c, i))
    if popl == 'none' and ener != 'none' and evtt == 'none' and back == 'none':
        for i in indxenertemp:
            liststrgvarb.append(strgvarb + 'en%02d' % i)
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none':
        liststrgvarb.append(strgvarb)
    
    return liststrgvarb


def setp_varbvalu(gdat, strgvarb, valu, popl='none', ener='none', evtt='none', back='none', isfr='none', strgmodl=None):
    
    if strgmodl is None:
        if gdat.datatype == 'mock':
            liststrgmodl = ['true', 'fitt']
        else:
            liststrgmodl = ['fitt']
    else:
        liststrgmodl = [strgmodl]
    liststrgmodl = deepcopy(liststrgmodl)
    for strgmodltemp in liststrgmodl:
        liststrgvarb = setp_namevarbsing(gdat, strgmodltemp, strgvarb, popl, ener, evtt, back, isfr)
        for strgvarbtemp in liststrgvarb:
            setp_varbiter(gdat, strgmodltemp, strgvarbtemp, valu)


def setp_varbiter(gdat, strgmodltemp, strgvarbtemp, valu):
    
    try:
        valutemp = getattr(gdat, strgvarbtemp)
        if valutemp is None:
            raise
        setattr(gdat, strgmodltemp + strgvarbtemp, valutemp)
    except:
        try:
            valutemp = getattr(gdat, strgmodltemp + strgvarbtemp)
            if valutemp is None:
                raise
        except:
            setattr(gdat, strgmodltemp + strgvarbtemp, valu)
    

def setp_varblimt(gdat, strgvarb, listvalu, typelimt='minmmaxm', popl='none', ener='none', evtt='none', back='none', isfr='none', strgmodl=None):
    
    if strgmodl is None:
        if gdat.datatype == 'mock':
            liststrgmodl = ['true', 'fitt']
        else:
            liststrgmodl = ['fitt']
    else:
        liststrgmodl = [strgmodl]
    
    for strgmodltemp in liststrgmodl:

        liststrgvarb = setp_namevarbsing(gdat, strgmodltemp, strgvarb, popl, ener, evtt, back, isfr)

        for strgvarbtemp in liststrgvarb:

            if typelimt == 'minmmaxm':
                setp_varbiter(gdat, strgmodltemp, 'minm' + strgvarbtemp, listvalu[0])
                setp_varbiter(gdat, strgmodltemp, 'maxm' + strgvarbtemp, listvalu[1])
            else:
                setp_varbiter(gdat, strgmodltemp, 'mean' + strgvarbtemp, listvalu[0])
                setp_varbiter(gdat, strgmodltemp, 'stdv' + strgvarbtemp, listvalu[1])
                
                # set np.minimum and np.maximum for Gaussian distributed variables
                meanpara = getattr(gdat, strgmodltemp + 'mean' + strgvarbtemp)
                stdp = getattr(gdat, strgmodltemp + 'stdv' + strgvarbtemp)
                setattr(gdat, strgmodltemp + 'minm' + strgvarbtemp, meanpara - stdp * 5)
                setattr(gdat, strgmodltemp + 'maxm' + strgvarbtemp, meanpara + stdp * 5)
 

def intp_sinc(gdat, lgal, bgal):

    intpsinc = 4. * gdat.numbsidepsfn**2 * np.sum(gdat.temppsfn * sinc(gdat.numbsidepsfn * (gdat.gridpsfnlgal + lgal) - gdat.gridpsfnlgal) * \
                                                                                                      sinc(gdat.numbsidepsfn * (gdat.gridpsfnbgal + bgal) - gdat.gridpsfnbgal))

    return intpsinc


def retr_fluxbrgt(gdat, lgal, bgal, flux):

    if lgal.size == 0:
        fluxbrgt = np.array([0.])
        fluxbrgtassc = np.array([0.])
    else:
        indxbrgt = np.argmax(flux)
        fluxbrgt = flux[indxbrgt]

    return fluxbrgt, fluxbrgtassc


def init_figr(gdat, gdatmodi, strgpdfn, strgplot, strgstat, strgmodl, indxenerplot, indxevttplot, indxpoplplot, intreval=False):
    
    if intreval:
        figrsize = [10, 10]
    else:
        figrsize = (gdat.sizeimag, gdat.sizeimag)
    figr, axis = plt.subplots(figsize=figrsize)
    
    if intreval:
        axis.set_position([0.5, 0.1, 0.6, 0.8])
        path = ''
    else:
        
        if indxenerplot is None:
            strgener = ''
        else:
            strgener = 'en%02d' % gdat.indxenerincl[indxenerplot]
        
        if indxevttplot is None:
            strgevtt = ''
        elif indxevttplot == -1:
            strgevtt = 'evtA'
        else:
            strgevtt = 'evt%d' % gdat.indxevttincl[indxevttplot]
        
        if indxpoplplot == -1:
            strgpopl = 'popA'
        else:
            strgpopl = 'pop%d' % indxpoplplot

        nameplot = '%s%s%s%s' % (strgplot, strgener, strgevtt, strgpopl)
   
        path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, nameplot)
    
    axis.set_xlabel(gdat.labllgaltotl)
    axis.set_ylabel(gdat.lablbgaltotl)
    titl = ''
    if indxenerplot is not None and gdat.numbener > 1 and strgplot.endswith('cnts'):
        titl = gdat.strgener[indxenerplot]
    if indxevttplot is not None and gdat.numbevtt > 1 and strgplot.endswith('cnts'):
        titl += ' ' + gdat.strgevtt[indxevttplot]
    axis.set_title(titl)

    return figr, axis, path


def draw_frambndr(gdat, axis):
    
    outr = max(gdat.frambndrmodl, gdat.frambndrdata)
    axis.set_xlim([-outr, outr])
    axis.set_ylim([-outr, outr])
    innr = min(gdat.frambndrmodl, gdat.frambndrdata)
    axis.axvline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axvline(-innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(-innr, ls='--', alpha=gdat.alphbndr, color='black')


def retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot=None, indxevttplot=-1, tdim=False, imag=None):
    
    vmin = getattr(gdat, 'minmscal' + strgcbar)
    vmax = getattr(gdat, 'maxmscal' + strgcbar)
    cmap = getattr(gdat, 'cmap' + strgcbar) 
    scal = getattr(gdat, 'scal' + strgcbar) 

    draw_frambndr(gdat, axis)
    
    # take the relevant energy and PSF bins
    if indxenerplot is not None:
        if indxevttplot == -1:
            maps = np.sum(maps[indxenerplot, ...], axis=1)
        else:
            maps = maps[indxenerplot, :, indxevttplot]
    
    # project the map to 2D
    if gdat.pixltype == 'heal':
        maps = tdpy.util.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                    minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                    minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata)
    
    if gdat.pixltype == 'cart':
        shap = [gdat.numbsidecart] + list(maps.shape)
        shap[1] = gdat.numbsidecart
        shapflat = list(maps.shape)
        shapflat[0] = gdat.numbpixlfull
        mapstemp = np.zeros(shapflat)
        if maps.size == gdat.indxpixlrofi.size:
            mapstemp[gdat.indxpixlrofi, ...] = maps
        else:
            mapstemp[:, ...] = maps
        maps = mapstemp.reshape(shap).swapaxes(0, 1)

    # temp -- this is needed to bring the Fermi-LAT map to the right direction
    #maps = fliplr(maps)

    # rescale the map
    if scal == 'asnh':
        maps = np.arcsinh(maps)
    if scal == 'logt':
        maps = np.log10(maps)
    if imag is None:
        imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest', vmin=vmin, vmax=vmax, alpha=gdat.alphmaps)
        return imag
    else:
        imag.set_data(maps)
    

def make_cbar(gdat, axis, imag, indxenerplot=None, tick=None, labl=None):

    # make a color bar
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05, aspect=15)
    if tick is not None and labl is not None:
        cbar.set_ticks(tick)
        cbar.set_ticklabels(labl)
    
    return cbar


def make_legdmaps(gdat, strgstat, strgmodl, axis, mosa=False, assc=False):
    
    # transdimensional elements
    if strgmodl == 'fitt' and (strgstat == 'pdfn' and gdat.boolcondcatl or strgstat == 'this') and gdat.fittnumbtrap > 0:
        for l in gdat.fittindxpopl:
            colr = retr_colr(gdat, strgstat, strgmodl, l)
            if strgstat == 'pdfn':
                labl = 'Condensed %s %s' % (gdat.fittlegd, gdat.fittlegdpopl[l])
            else:
                labl = 'Sample %s %s' % (gdat.fittlegd, gdat.fittlegdpopl[l])
            if not gdat.fittmaxmnumbelempopl[l] == 0:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                                                        label=labl, marker=gdat.fittlistelemmrkr[l], lw=gdat.mrkrlinewdth, color=colr)
    
    for q in gdat.indxrefr:
        if not np.amax(gdat.refrnumbelem[q]) == 0:
            if assc:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                       label=gdat.refrlegdhits[q], marker=gdat.refrlistelemmrkrhits[q], lw=gdat.mrkrlinewdth, color=gdat.refrcolrelem[q])
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                       label=gdat.refrlegdmiss[q], marker=gdat.refrlistelemmrkrmiss[q], lw=gdat.mrkrlinewdth, color=gdat.refrcolrelem[q])
            else:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                       label=gdat.refrlegdelem[q], marker=gdat.refrlistelemmrkrmiss[q], lw=gdat.mrkrlinewdth, color=gdat.refrcolrelem[q])
    
    # fixed-dimensional objects
    if strgmodl == 'fitt':
        if gdat.fittlensmodltype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                           label='%s Source' % gdat.fittlegd, marker='<', lw=gdat.mrkrlinewdth, color=gdat.fittcolr)
        
        if gdat.fitthostemistype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                           label='%s Host' % gdat.fittlegd, marker='s', lw=gdat.mrkrlinewdth, color=gdat.fittcolr)
    
    if gdat.datatype == 'mock':
        if gdat.truelensmodltype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                label='%s Source' % gdat.refrlegd, marker='>', lw=gdat.mrkrlinewdth, color=gdat.refrcolr)
        
        if gdat.truehostemistype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                label='%s Host' % gdat.refrlegd, marker='D', lw=gdat.mrkrlinewdth, color=gdat.refrcolr)
    
    temphand, temp = axis.get_legend_handles_labels()
    numblabl = len(temp)
    
    if numblabl == 4:
        numbcols = 2
    else:
        numbcols = 3
    if mosa:
        axis.legend(bbox_to_anchor=[1., 1.15], loc='center', ncol=numbcols)
    else:
        axis.legend(bbox_to_anchor=[0.5, 1.15], loc='center', ncol=numbcols)
        

def supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxpoplplot=-1, assc=False):
    
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    
    if lensmodltype != 'none' or hostemistype != 'none':
        indxsersfgrd = getattr(gdat, strgmodl + 'indxsersfgrd')
    
    # associations with the reference elements
    for q in gdat.indxrefr:
        if gdat.refrnumbelem[q] > 0:
            if indxpoplplot == -1:
                listindxpoplplot = gdat.fittindxpopl
            else:
                listindxpoplplot = [indxpoplplot]
            for l in listindxpoplplot:
                # temp -- backcomp
                try:
                    reframpl = getattr(gdat, 'refr' + gdat.refrlistnamefeatampl[q])
                except:
                    reframpl = getattr(gdat, 'refr' + gdat.listnamefeatamplrefr[q])
                if reframpl is None:
                    mrkrsize = full(gdat.refrnumbelem[q], 5.)
                else:
                    # temp -- backcomp
                    try:
                        mrkrsize = retr_mrkrsize(gdat, reframpl[q][0, :], gdat.refrlistnamefeatampl[q])
                    except:
                        mrkrsize = retr_mrkrsize(gdat, reframpl[q][0, :], gdat.listnamefeatamplrefr[q])
                lgal = np.copy(gdat.refrlgal[q][0, :])
                bgal = np.copy(gdat.refrbgal[q][0, :])
                numbelem = int(gdat.refrnumbelem[q])
                
                if gdatmodi is not None and numbtrap > 0 and assc:   
                    ### hit
                    indx = gdatmodi.thisindxelemrefrasschits[q][l]
                    if indx.size > 0:
                        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, label=gdat.refrlegdhits, \
                                                                      marker=gdat.refrlistelemmrkrhits[q], lw=gdat.mrkrlinewdth, color=gdat.refrcolrelem[q])
                    ### missed
                    indx = gdatmodi.thisindxelemrefrasscmiss[q][l]
                else:
                    indx = np.arange(lgal.size)
                if indx.size > 0: 
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, facecolor='none', \
                                                             label=gdat.refrlegdmiss, marker=gdat.refrlistelemmrkrmiss[q], lw=gdat.mrkrlinewdth, color=gdat.refrcolrelem[q])
        
            sizexoff = gdat.maxmgangdata * 0.05 * gdat.anglfact
            sizeyoff = gdat.maxmgangdata * 0.05 * gdat.anglfact
            if 'etag' in gdat.refrliststrgfeat[q]:
                for k in range(indx.size):
                    axis.text(gdat.anglfact * lgal[indx[k]] + sizexoff, gdat.anglfact * bgal[indx[k]] + sizeyoff, gdat.refretag[q][indx[k]], \
                                                                                            verticalalignment='center', horizontalalignment='center', \
                                                                                                             color='red', fontsize=1)

    # temp -- generalize this to input refrlgalhost vs.
    if gdat.datatype == 'mock':
        ## host galaxy position
        if hostemistype != 'none':
            for e in indxsersfgrd:
                lgalhost = gdat.truesamp[getattr(gdat, 'trueindxfixplgalhostisf%d' % (e))]
                bgalhost = gdat.truesamp[getattr(gdat, 'trueindxfixpbgalhostisf%d' % (e))]
                axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, facecolor='none', alpha=0.7, \
                                             label='%s Host %d' % (gdat.refrlegd, e), s=300, marker='D', lw=gdat.mrkrlinewdth, color=gdat.refrcolr)
        if lensmodltype != 'none':
            ## host galaxy Einstein radius
            for e in gdat.trueindxsersfgrd:
                truelgalhost = gdat.truesamp[getattr(gdat, 'trueindxfixplgalhostisf%d' % (e))]
                truebgalhost = gdat.truesamp[getattr(gdat, 'trueindxfixpbgalhostisf%d' % (e))]
                truebeinhost = gdat.truesamp[getattr(gdat, 'trueindxfixpbeinhostisf%d' % (e))]
                axis.add_patch(plt.Circle((gdat.anglfact * truelgalhost, \
                                           gdat.anglfact * truebgalhost), \
                                           gdat.anglfact * truebeinhost, \
                                           edgecolor=gdat.refrcolr, facecolor='none', lw=gdat.mrkrlinewdth))
            
        if lensmodltype != 'none':
            ## source galaxy position
            axis.scatter(gdat.anglfact * gdat.truesamp[gdat.trueindxfixplgalsour], \
                                                    gdat.anglfact * gdat.truesamp[gdat.trueindxfixpbgalsour], \
                                                    facecolor='none', \
                                                    alpha=0.7, \
                                                    #alpha=gdat.alphelem, \
                                                    label='%s Source' % gdat.refrlegd, s=300, marker='>', lw=gdat.mrkrlinewdth, color=gdat.refrcolr)
        
    # model catalog
    if indxpoplplot == -1:
        listindxpoplplot = gdat.fittindxpopl
    else:
        listindxpoplplot = [indxpoplplot]
    for l in listindxpoplplot:
        if gdatmodi is not None:
            if gdat.fittnumbtrap > 0:
                colr = retr_colr(gdat, strgstat, strgmodl, l)
                mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissamp[gdatmodi.thisindxsampcomp[gdat.fittnamefeatampl[l]][l]], gdat.fittnamefeatampl[l])
                if 'lgal' in gdatmodi.thisindxsampcomp:
                    lgal = gdatmodi.thissamp[gdatmodi.thisindxsampcomp['lgal'][l]]
                    bgal = gdatmodi.thissamp[gdatmodi.thisindxsampcomp['bgal'][l]]
                else:
                    gang = gdatmodi.thissamp[gdatmodi.thisindxsampcomp['gang'][l]]
                    aang = gdatmodi.thissamp[gdatmodi.thisindxsampcomp['aang'][l]]
                    lgal, bgal = retr_lgalbgal(gang, aang)
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphelem, label='Sample', marker=gdat.fittlistelemmrkr[l], \
                                                                                                                           lw=gdat.mrkrlinewdth, color=colr)

            ## source
            if lensmodltype != 'none':
                lgalsour = gdatmodi.thissamp[gdat.fittindxfixplgalsour]
                bgalsour = gdatmodi.thissamp[gdat.fittindxfixpbgalsour]
                axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, facecolor='none', \
                                                  alpha=gdat.alphelem, \
                                                  label='%s Source' % gdat.fittlegd, s=300, marker='<', lw=gdat.mrkrlinewdth, color=gdat.fittcolr)
    
            if hostemistype != 'none':
                ## host
                lgalhost = [[] for e in indxsersfgrd]
                bgalhost = [[] for e in indxsersfgrd]
                for e in indxsersfgrd:
                    lgalhost[e] = gdatmodi.thissamp[getattr(gdat, strgmodl + 'indxfixplgalhostisf%d' % (e))]
                    bgalhost[e] = gdatmodi.thissamp[getattr(gdat, strgmodl + 'indxfixpbgalhostisf%d' % (e))]
                    axis.scatter(gdat.anglfact * lgalhost[e], gdat.anglfact * bgalhost[e], facecolor='none', \
                                                     alpha=gdat.alphelem, \
                                                     label='%s Host' % gdat.fittlegd, s=300, marker='s', lw=gdat.mrkrlinewdth, color=gdat.fittcolr)
                    if lensmodltype != 'none':
                        beinhost = gdatmodi.thissamp[getattr(gdat, strgmodl + 'indxfixpbeinhostisf%d' % (e))]
                        try:
                            axis.add_patch(plt.Circle((gdat.anglfact * lgalhost[e], gdat.anglfact * bgalhost[e]), \
                                                       gdat.anglfact * beinhost, edgecolor=gdat.fittcolr, facecolor='none', \
                                                       lw=gdat.mrkrlinewdth, ls='--'))
                        except:
                            print 'Circle failed...'
                            print 'gdat.anglfact'
                            print gdat.anglfact
                            print 'beinhost'
                            print beinhost
                            print 'lgalhost'
                            print lgalhost
                            print 'bgalhost'
                            print bgalhost
                            print 'gdat.fittcolr'
                            print gdat.fittcolr
                            print 'gdat.mrkrlinewdth'
                            print gdat.mrkrlinewdth
                            pass
                            raise Exception('')
                
    # temp
    if strgstat == 'pdfn' and gdat.boolcondcatl and gdat.fittnumbtrap > 0:
        lgal = np.zeros(gdat.numbprvlhigh)
        bgal = np.zeros(gdat.numbprvlhigh)
        ampl = np.zeros(gdat.numbprvlhigh)
        cntr = 0
        for r in gdat.indxstkscond:
            if r in gdat.indxprvlhigh:
                lgal[cntr] = gdat.dictglob['poststkscond'][r]['lgal'][0]
                bgal[cntr] = gdat.dictglob['poststkscond'][r]['bgal'][0]
                # temp -- this does not allow sources with different spectra to be assigned to the same stacked sample
                ampl[cntr] = gdat.dictglob['poststkscond'][r][gdat.fittnamefeatampl[l]][0]
                cntr += 1
        mrkrsize = retr_mrkrsize(gdat, ampl, gdat.fittnamefeatampl[l])
        
        colr = retr_colr(gdat, strgstat, strgmodl, l)
        axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, \
                                    label='Condensed', marker=gdat.fittlistelemmrkr[l], color='black', lw=gdat.mrkrlinewdth)
        for r in gdat.indxstkscond:
            lgal = np.array([gdat.dictglob['liststkscond'][r]['lgal']])
            bgal = np.array([gdat.dictglob['liststkscond'][r]['bgal']])
            axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, \
                                                marker=gdat.fittlistelemmrkr[l], color='black', alpha=0.1, lw=gdat.mrkrlinewdth)


def retr_colr(gdat, strgstat, strgmodl, indxpopl=None):
    
    if strgmodl == 'true':
        if indxpopl is None:
            colr = gdat.refrcolr
        else:
            colr = gdat.refrcolrelem[indxpopl]
    if strgmodl == 'fitt':
        if strgstat == 'this' or strgstat == 'pdfn':
            if indxpopl is None:
                colr = gdat.fittcolr
            else:
                colr = gdat.fittcolrelem[indxpopl]
        if strgstat == 'mlik':
            colr = 'r'
    
    return colr


def retr_levipost(listllik):
    
    minmlistllik = np.amin(listllik)
    levipost = np.log(np.mean(1. / np.exp(listllik - minmlistllik))) + minmlistllik
    
    return levipost


def retr_infofromlevi(pmeallik, levi):
    
    info = pmeallik - levi

    return info


def retr_jcbn():
    
    fluxpare, lgalpare, bgalpare, fluxauxi, lgalauxi, bgalauxi = sympy.symbols('fluxpare lgalpare bgalpare fluxauxi lgalauxi bgalauxi')
    
    matr = sympy.Matrix([[ fluxpare,      fluxauxi, 0,            0, 0,            0], \
                         [-fluxpare, 1  - fluxauxi, 0,            0, 0,            0], \
                         [-lgalauxi,             0, 1, 1 - fluxauxi, 0,            0], \
                         [-lgalauxi,             0, 1,    -fluxauxi, 0,            0], \
                         [-bgalauxi,             0, 0,            0, 1, 1 - fluxauxi], \
                         [-bgalauxi,             0, 0,            0, 1,    -fluxauxi]])

    jcbn = matr.det()
    print jcbn

    return jcbn

# f1 = uf f0
# f2 = (1 - uf) f0
# x1 = x0 + (1 - uf) ux
# x2 = x0 - uf ux
# y1 = y0 + (1 - uf) uy
# y2 = y0 - uf uy

# f1/uf f1/f0 f1/x0 f1/ux f1/y0 f1/uy
# f2/uf f2/f0 f2/x0 f2/ux f2/y0 f2/uy
# x1/uf x1/f0 x1/x0 x1/ux x1/y0 x1/uy
# x2/uf x2/f0 x2/x0 x2/ux x2/y0 x2/uy
# y1/uf y1/f0 y1/x0 y1/ux y1/y0 y1/uy
# y2/uf y2/f0 y2/x0 y2/ux y2/y0 y2/uy

#  f0     uf 0      0 0      0
# -f0 1 - uf 0      0 0      0
# -ux      0 1 1 - uf 0      0
# -ux      0 1    -uf 0      0
# -uy      0 0      0 1 1 - uf
# -uy      0 0      0 1    -uf

# f0
#retr_jcbn()

def retr_angldist(gdat, lgalfrst, bgalfrst, lgalseco, bgalseco):
    
    # temp -- heal does not work when the dimension of lgalfrst is 1
    if gdat.pixltype == 'heal':
        dir1 = np.array([lgalfrst, bgalfrst])
        dir2 = np.array([lgalseco, bgalseco])
        angldist = hp.rotator.angdist(dir1, dir2)
    else:
        angldist = np.sqrt((lgalfrst - lgalseco)**2 + (bgalfrst - bgalseco)**2)

    return angldist


def retr_deflextr(gdat, indxpixlelem, sher, sang):
    
    factcosi = sher * np.cos(2. * sang)
    factsine = sher * np.cos(2. * sang)
    defllgal = factcosi * gdat.lgalgrid[indxpixlelem] + factsine * gdat.bgalgrid[indxpixlelem]
    deflbgal = factsine * gdat.lgalgrid[indxpixlelem] - factcosi * gdat.bgalgrid[indxpixlelem]
    
    return np.vstack((defllgal, deflbgal)).T 


def readfile(path):

    print 'Reading %s...' % path

    filepick = open(path + '.p', 'rb')
    filearry = h5py.File(path + '.h5', 'r')
    gdattemptemp = cPickle.load(filepick)
    
    for attr in filearry:
        setattr(gdattemptemp, attr, filearry[attr][()])

    filepick.close()
    filearry.close()
    
    if 'gdatfinl' in path or 'gdatinit' in path:
        if hasattr(gdattemptemp, 'edis') and gdattemptemp.edis is not None and hasattr(gdattemptemp, 'binsener'):
            gdattemptemp.edisintp = interp1d_pick(gdattemptemp.binsener, gdattemptemp.edis)
        gdattemptemp.adisobjt = interp1d_pick(gdattemptemp.redsintp, gdattemptemp.adisintp)
        gdattemptemp.redsfromdlosobjt = interp1d_pick(gdattemptemp.adisintp * gdattemptemp.redsintp, gdattemptemp.redsintp)
    
    return gdattemptemp


def writfile(gdattemp, path):
    
    filepick = open(path + '.p', 'wb')
    filearry = h5py.File(path + '.h5', 'w')
    
    gdattemptemp = tdpy.util.gdatstrt()
    for attr, valu in gdattemp.__dict__.iteritems():
        if attr.endswith('psfnintp'):
            continue

        if isinstance(valu, np.ndarray) and valu.dtype != np.dtype('O') or isinstance(valu, str) or \
                                       isinstance(valu, float) or isinstance(valu, bool) or isinstance(valu, int) or isinstance(valu, np.float):
            filearry.create_dataset(attr, data=valu)
        else:
            # temp -- make sure interpolation objects are not written.
            if attr != 'adisobjt' and attr != 'redsfromdlosobjt' and attr != 'edisintp':
                setattr(gdattemptemp, attr, valu)
    
    print 'Writing to %s...' % path

    cPickle.dump(gdattemptemp, filepick, protocol=cPickle.HIGHEST_PROTOCOL)
    filepick.close()
    filearry.close()
   

def retr_deflcutf(angl, defs, asca, acut, asym=False):

    fracanglasca = angl / asca
    
    deflcutf = defs / fracanglasca
    
    # second term in the NFW deflection profile
    fact = np.ones_like(fracanglasca)
    indxlowr = np.where(fracanglasca < 1.)[0]
    indxuppr = np.where(fracanglasca > 1.)[0]
    fact[indxlowr] = arccosh(1. / fracanglasca[indxlowr]) / np.sqrt(1. - fracanglasca[indxlowr]**2)
    fact[indxuppr] = np.arccos(1. / fracanglasca[indxuppr]) / np.sqrt(fracanglasca[indxuppr]**2 - 1.)
    
    if asym:
        deflcutf *= np.log(fracanglasca / 2.) + fact
    else:
        fracacutasca = acut / asca
        factcutf = fracacutasca**2 / (fracacutasca**2 + 1)**2 * ((fracacutasca**2 + 1. + 2. * (fracanglasca**2 - 1.)) * fact + \
                np.pi * fracacutasca + (fracacutasca**2 - 1.) * np.log(fracacutasca) + np.sqrt(fracanglasca**2 + fracacutasca**2) * (-pi + (fracacutasca**2 - 1.) / fracacutasca * \
                np.log(fracanglasca / (np.sqrt(fracanglasca**2 + fracacutasca**2) + fracacutasca))))
        deflcutf *= factcutf
       
    return deflcutf


def initchro(gdat, gdatmodi, name):
    
    if gdatmodi is not None:    
        setattr(gdatmodi, 'thischro' + name, gdat.functime())
    

def stopchro(gdat, gdatmodi, name):
    
    if gdatmodi is not None:    
        setattr(gdatmodi, 'thischro' + name, gdat.functime() - getattr(gdatmodi, 'thischro' + name))


def retr_defl(gdat, indxpixlelem, lgal, bgal, angllens, ellp=None, angl=None, rcor=None, asca=None, acut=None):
    
    # translate the grid
    lgaltran = gdat.lgalgrid[indxpixlelem] - lgal
    bgaltran = gdat.bgalgrid[indxpixlelem] - bgal
    
    if acut is not None:
        defs = angllens
        angl = np.sqrt(lgaltran**2 + bgaltran**2)
        defl = retr_deflcutf(angl, defs, asca, acut)
        defllgal = lgaltran / angl * defl
        deflbgal = bgaltran / angl * defl

    else:
        bein = angllens

        # rotate the grid
        lgalrttr = np.cos(angl) * lgaltran - np.sin(angl) * bgaltran
        bgalrttr = np.sin(angl) * lgaltran + np.cos(angl) * bgaltran
        
        axisrati = 1. - ellp
        facteccc = np.sqrt(1. - axisrati**2)
        factrcor = np.sqrt(axisrati**2 * lgalrttr**2 + bgalrttr**2)
        defllgalrttr = bein * axisrati / facteccc *  np.arctan(facteccc * lgalrttr / factrcor)
        deflbgalrttr = bein * axisrati / facteccc * np.arctanh(facteccc * bgalrttr / factrcor)
        
        # totate back vector to original basis
        defllgal = np.cos(angl) * defllgalrttr + np.sin(angl) * deflbgalrttr
        deflbgal = -np.sin(angl) * defllgalrttr + np.cos(angl) * deflbgalrttr
   
    defl = np.vstack((defllgal, deflbgal)).T
    
    return defl


def retr_lpriselfdist(gdat, strgmodl, feat, strgfeat):
    
    minm = getattr(gdat, 'minm' + strgfeat)
    maxm = getattr(gdat, 'maxm' + strgfeat)
    
    lpri = np.sum(np.log(pdfn_self(feat, minm, maxm)))
    
    return lpri


def retr_lprilogtdist(gdat, strgmodl, feat, strgfeat):
    
    minm = getattr(gdat, 'minm' + strgfeat)
    maxm = getattr(gdat, 'maxm' + strgfeat)
    
    lpri = np.sum(np.log(pdfn_logt(feat, minm, maxm)))
    
    return lpri


def retr_lpripowrdist(gdat, strgmodl, feat, strgfeat, samp, l):
    
    minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
    maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
    distslop = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
    
    lpri = np.sum(np.log(pdfn_powr(feat, minm, maxm, distslop)))
    
    return lpri


def retr_lpridpowdist(gdat, strgmodl, feat, strgfeat, samp, l):
    
    minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
    maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
    brek = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distbrek')[l]]
    distsloplowr = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distsloplowr')[l]]
    distslopuppr = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslopuppr')[l]]
    
    lpri = np.sum(np.log(pdfn_dpow(feat, minm, maxm, brek, distsloplowr, distslopuppr)))
    
    return lpri


def retr_lprigausdist(gdat, strgmodl, feat, strgfeat, samp, l):
    
    distmean = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distmean')[l]]
    diststdv = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'diststdv')[l]]
    
    lpri = np.sum(np.log(pdfn_gaus(feat, distmean, diststdv)))
    
    return lpri


def retr_lpriigamdist(gdat, strgmodl, feat, strgfeat, samp, l):
    
    distslop = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
    cutf = getattr(gdat, strgmodl + 'cutf' + strgfeat)
    
    lpri = np.sum(np.log(pdfn_igam(feat, distslop, cutf)))

    return lpri


def traptdim(gdat, arry):
    
    s1 = arry[0, 0] + arry[-1, 0] + arry[0, -1] + arry[-1, -1]
    s2 = np.sum(arry[1:-1, 0]) + np.sum(arry[1:-1, -1]) + np.sum(arry[0, 1:-1]) + np.sum(arry[-1, 1:-1])
    s3 = np.sum(arry[1:-1, 1:-1])
    summ = (s1 + 2*s2 + 4*s3) * gdat.apix
    
    return summ


def retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons=None):
    
    pdfnspatprio = pdfnspatpriotemp
    if spatdistcons is not None:
        pdfnspatprio += spatdistcons

    summ = traptdim(gdat, pdfnspatprio)
    pdfnspatprio /= summ
    lpdfspatprio = np.log(pdfnspatprio)
    lpdfspatprioobjt = sp.interpolate.RectBivariateSpline(gdat.binsbgalcart, gdat.binslgalcart, lpdfspatprio)
    
    return lpdfspatprio, lpdfspatprioobjt


def retr_strgpfix(strgstat, strgmodl, strgpdfn='post', strgmome='pmea'):
    
    if strgstat == 'pdfn':
        strgpfix = strgmome + strgpdfn
    else:
        if strgmodl == 'true':
            if strgstat == 'this':
                strgpfix = 'true'
            else:
                strgpfix = 'nexttrue'
        else:
            strgpfix = strgstat
    
    return strgpfix

        
def retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl):
    
    if strgmodl == 'true' or strgmodl == 'fitt' and (strgstat == 'mlik' or strgstat == 'pdfn'):
        gdatobjt = gdat
    else:
        gdatobjt = gdatmodi

    return gdatobjt


def proc_samp(gdat, gdatmodi, strgstat, strgmodl, raww=False, fast=False, boolinit=False):
   
    initchro(gdat, gdatmodi, 'pars')
    if gdat.verbtype > 1:
        print 'proc_samp()'
        print 'strgstat'
        print strgstat
    
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl)
    fittnumbelemzero = getattr(gdat, strgmodl + 'numbelemzero')
    # temp
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    boolbfun = getattr(gdat, strgmodl + 'boolbfun')
    if numbtrap > 0:
        indxfixpmeanelem = getattr(gdat, strgmodl + 'indxfixpmeanelem')
        boolelemlens = getattr(gdat, strgmodl + 'boolelemlens')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        elemtype = getattr(gdat, strgmodl + 'elemtype')
        boolelemlghtanyy = getattr(gdat, strgmodl + 'boolelemlghtanyy')
        boolelemsbrtdfnc = getattr(gdat, strgmodl + 'boolelemsbrtdfnc')
        spatdisttype = getattr(gdat, strgmodl + 'spatdisttype')
        boolelemlght = getattr(gdat, strgmodl + 'boolelemlght')
        boolelemdeflsubh = getattr(gdat, strgmodl + 'boolelemdeflsubh')
        boolelemdeflsubhanyy = getattr(gdat, strgmodl + 'boolelemdeflsubhanyy')
        boolelemsbrtdfncanyy = getattr(gdat, strgmodl + 'boolelemsbrtdfncanyy')
        boolelemsbrtextsbgrd = getattr(gdat, strgmodl + 'boolelemsbrtextsbgrd')
        indxpopllens = getattr(gdat, strgmodl + 'indxpopllens')
        indxcompampl = getattr(gdat, strgmodl + 'indxcompampl')
        boolelemsbrtextsbgrdanyy = getattr(gdat, strgmodl + 'boolelemsbrtextsbgrdanyy')
        boolelempsfn = getattr(gdat, strgmodl + 'boolelempsfn')
        boolelempsfnanyy = getattr(gdat, strgmodl + 'boolelempsfnanyy')
        boolelemspat = getattr(gdat, strgmodl + 'boolelemspat')
        boolelemlghtspat = getattr(gdat, strgmodl + 'boolelemlghtspat')
        if boolelemlghtanyy:
            minmflux = getattr(gdat, strgmodl + 'minmflux')
        if 'lens' in elemtype:
            minmdefs = getattr(gdat, strgmodl + 'minmdefs')
        if 'clus' in elemtype or 'clusvari' in elemtype:
            minmnobj = getattr(gdat, strgmodl + 'minmnobj')
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        numbcomp = getattr(gdat, strgmodl + 'numbcomp')
        spectype = getattr(gdat, strgmodl + 'spectype')
        liststrgfeatdefa = getattr(gdat, strgmodl + 'liststrgfeatdefa')
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        liststrgfeatodimtotl = getattr(gdat, strgmodl + 'liststrgfeatodimtotl')
        indxgridpopl = getattr(gdat, strgmodl + 'indxgridpopl')
        liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
        indxpoplgrid = getattr(gdat, strgmodl + 'indxpoplgrid')
        listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
        liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
        maxmnumbcomp = getattr(gdat, strgmodl + 'maxmnumbcomp')
        liststrgfeateval = getattr(gdat, strgmodl + 'liststrgfeateval')
        elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
    
        if gdat.verbtype > 2:
            print 'spectype'
            print spectype
    
    backtype = getattr(gdat, strgmodl + 'backtype')
    psfntype = getattr(gdat, strgmodl + 'psfntype')
    numblpri = getattr(gdat, strgmodl + 'numblpri')
    #convdiff = getattr(gdat, strgmodl + 'convdiff')
    convdiffanyy = getattr(gdat, strgmodl + 'convdiffanyy')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    listnamediff = getattr(gdat, strgmodl + 'listnamediff')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')
    
    strgpfix = retr_strgpfix(strgstat, strgmodl)
    strgpfixthis = retr_strgpfix('this', strgmodl)
    strgpfixnext = retr_strgpfix('next', strgmodl)

    # grab the sample vector
    samp = getattr(gdatobjt, strgpfix + 'samp')

    indxpara = np.arange(samp.size) 

    if gdat.booldiagmode:
        if not np.isfinite(samp).all():
            print 'samp'
            indx = np.where(np.logical_not(np.isfinite(samp)))
            print 'samp[indx]'
            print samp[indx]
            namepara = getattr(gdat, strgmodl + 'namepara')
            print namepara[indx]
            raise Exception('')

    if psfnevaltype != 'none':
        psfp = samp[getattr(gdat, strgmodl + 'indxfixppsfp')]

    bacp = samp[getattr(gdat, strgmodl + 'indxfixpbacp')]
    
    if numbtrap > 0:
        indxelemfull = list(getattr(gdatobjt, strgpfix + 'indxelemfull'))
        if gdat.verbtype > 2:
            print 'indxelemfull'
            print indxelemfull
        # temp -- this may slow down execution
        indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, strgmodl)

        setattr(gdatobjt, strgpfix + 'indxsampcomp', indxsampcomp)
        indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
    
        numbelempopl = np.empty(numbpopl)
        numbelem = [[] for l in indxpopl]
        for l in indxpopl:
            numbelem[l] = samp[indxfixpnumbelem[l]].astype(int)
            numbelempopl[l] = np.sum(numbelem[l])
        numbelemtotl = np.sum(numbelempopl) 
       
        dictelem = [[] for l in indxpopl]
        for l in indxpopl:
            dictelem[l] = dict()
            for strgfeat in liststrgfeatdefa:
                dictelem[l][strgfeat] = []
            for strgcomp in liststrgcomp[l]:
                dictelem[l][strgcomp] = samp[indxsampcomp[strgcomp][l]]
                if gdat.booldiagmode:
                    if ((abs(samp[indxsampcomp[strgcomp][l]]) < 1e-100 ) & (abs(samp[indxsampcomp[strgcomp][l]]) > 0.)).any():
                        print 'l'
                        print l
                        print 'strgcomp'
                        print strgcomp
                        print 'samp[indxsampcomp[strgcomp][l]]'
                        summgene(samp[indxsampcomp[strgcomp][l]])
                        print 'Parameter is too small!'
                        raise Exception('')

                    if numbelem[l] != len(dictelem[l][strgcomp]):
                        print 'l'
                        print l
                        print 'numbelem[l]'
                        print numbelem[l]
                        print 'indxsampcomp'
                        print indxsampcomp
                        print 'strgcomp'
                        print strgcomp
                        print 'dictelem[l][strgcomp]'
                        summgene(dictelem[l][strgcomp])
                        raise Exception('')
    
        if gdat.booldiagmode:
            for l in indxpopl:
                for g, strgcomp in enumerate(liststrgcomp[l]):
                    if (listscalcomp[l][g] != 'gaus' and not listscalcomp[l][g].startswith('lnor')) and  \
                                        (listscalcomp[l][g] != 'expo' and (dictelem[l][strgcomp] < getattr(gdat, strgmodl + 'minm' + strgcomp)).any()) or \
                                        (dictelem[l][strgcomp] > getattr(gdat, strgmodl + 'maxm' + strgcomp)).any():
                        
                        print 'l'
                        print l
                        print 'strgcomp'
                        print strgcomp
                        print 'dictelem[l][strgcomp]'
                        print dictelem[l][strgcomp]
                        print 'getattr(gdat, strgmodl + minm + strgcomp)'
                        print getattr(gdat, strgmodl + 'minm' + strgcomp)
                        print 'getattr(gdat, strgmodl + maxm + strgcomp)'
                        print getattr(gdat, strgmodl + 'maxm' + strgcomp)
                        print '(dictelem[l][strgcomp] > getattr(gdat, strgmodl + maxm + strgcomp))'
                        print (dictelem[l][strgcomp] > getattr(gdat, strgmodl + 'maxm' + strgcomp))
                        print '(listscalcomp[l][g] != expo and (dictelem[l][strgcomp] < getattr(gdat, strgmodl + minm + strgcomp)).any())'
                        print (listscalcomp[l][g] != 'expo' and (dictelem[l][strgcomp] < getattr(gdat, strgmodl + 'minm' + strgcomp)).any())
                        print 'dictelem[l][strgcomp] < getattr(gdat, strgmodl + minm + strgcomp)'
                        print dictelem[l][strgcomp] < getattr(gdat, strgmodl + 'minm' + strgcomp)
                        raise Exception('')
                        print('Element parameter outside prior.')
                        #raise Exception('Element parameter outside prior.')
           
        # temp
        if gdat.booldiagmode:
            for l in indxpopl:
                if elemtype[l] == 'lens':
                    if gdat.variasca:
                        indx = np.where(samp[indxsampcomp['acut'][l]] < 0.)[0]
                        if indx.size > 0:
                            print 'Acut went negative'
                            samp[indxsampcomp['acut'][l]][indx] = 1e-3 * gdat.anglfact
                            raise Exception('')
                    if gdat.variacut:
                        indx = np.where(samp[indxsampcomp['asca'][l]] < 0.)[0]
                        if indx.size > 0:
                            print 'Asca went negative'
                            samp[indxsampcomp['asca'][l]][indx] = 1e-3 * gdat.anglfact
                            raise Exception('')
    
        for l in indxpopl:
            if elemtype[l].startswith('lght'):
                    
                # evaluate horizontal and vertical position for elements whose position is a power law in image-centric radius
                if spatdisttype[l] == 'glc3':
                    dictelem[l]['dlos'], dictelem[l]['lgal'], dictelem[l]['bgal'] = retr_glc3(dictelem[l]['dglc'], \
                                                                                                        dictelem[l]['thet'], dictelem[l]['phii'])
                
                if spatdisttype[l] == 'gangexpo':
                    dictelem[l]['lgal'], dictelem[l]['bgal'], = retr_lgalbgal(dictelem[l]['gang'], \
                                                                                                        dictelem[l]['aang'])
                    
                    if gdat.booldiagmode:
                        if numbelem[l] > 0:
                            if np.amin(dictelem[l]['lgal']) < getattr(gdat, strgmodl + 'minmlgal') or \
                               np.amax(dictelem[l]['lgal']) > getattr(gdat, strgmodl + 'maxmlgal') or \
                               np.amin(dictelem[l]['bgal']) < getattr(gdat, strgmodl + 'minmbgal') or \
                               np.amax(dictelem[l]['bgal']) > getattr(gdat, strgmodl + 'maxmbgal'):
                                raise Exception('Bad coordinates!')

                if spatdisttype[l] == 'los3':
                    dictelem[l]['dglc'], dictelem[l]['thet'], dictelem[l]['phii'] = retr_los3(dictelem[l]['dlos'], \
                                                                                                        dictelem[l]['lgal'], dictelem[l]['bgal'])

                # evaluate flux for pulsars
                if elemtype[l] == 'lghtpntspuls':
                    dictelem[l]['lumi'] = retr_lumipuls(dictelem[l]['geff'], dictelem[l]['magf'], dictelem[l]['per0'])
                if elemtype[l] == 'lghtpntsagnntrue':
                    dictelem[l]['reds'] = gdat.redsfromdlosobjt(dictelem[l]['dlos'])
                    dictelem[l]['lumi'] = dictelem[l]['lum0'] * (1. + dictelem[l]['reds'])**4
                if elemtype[l] == 'lghtpntspuls' or elemtype[l] == 'lghtpntsagnntrue':
                    dictelem[l]['flux'] = retr_flux(gdat, dictelem[l]['lumi'], dictelem[l]['dlos'])
                # evaluate spectra
                if elemtype[l].startswith('lghtline'):
                    if elemtype[l] == 'lghtlinevoig':
                        dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], elin=dictelem[l]['elin'], sigm=dictelem[l]['sigm'], \
                                                                                                          gamm=dictelem[l]['gamm'], spectype=spectype[l])
                    else:
                        dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], elin=dictelem[l]['elin'], edisintp=gdat.edisintp, spectype=spectype[l])
                else:
                    sindcolr = [dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                    dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], sind=dictelem[l]['sind'], curv=dictelem[l]['curv'], \
                                       expc=dictelem[l]['expc'], sindcolr=sindcolr, spectype=spectype[l])
        

    stopchro(gdat, gdatmodi, 'pars')
    
    ### loglikelihood
    initchro(gdat, gdatmodi, 'modl')
    
    if lensmodltype != 'none' or hostemistype != 'none':
        indxsersfgrdeval = getattr(gdat, strgmodl + 'indxsersfgrd')
    
    if lensmodltype != 'none':
        lgalsour = samp[getattr(gdat, strgmodl + 'indxfixplgalsour')]
        bgalsour = samp[getattr(gdat, strgmodl + 'indxfixpbgalsour')]
    
    if gdat.verbtype > 2:
        print 'Evaluating the likelihood...'
    
    # process a sample vector and the occupancy list to calculate secondary variables
    if lensmodltype != 'none':
        fluxsour = samp[getattr(gdat, strgmodl + 'indxfixpfluxsour')]
        if gdat.numbener > 1:
            sindsour = samp[getattr(gdat, strgmodl + 'indxfixpsindsour')]
        sizesour = samp[getattr(gdat, strgmodl + 'indxfixpsizesour')]
        ellpsour = samp[getattr(gdat, strgmodl + 'indxfixpellpsour')]
        anglsour = samp[getattr(gdat, strgmodl + 'indxfixpanglsour')]
    if lensmodltype != 'none' or hostemistype != 'none':
        indxsersfgrd = getattr(gdat, strgmodl + 'indxsersfgrd')
    if hostemistype != 'none':
        lgalhost = [[] for e in indxsersfgrd]
        bgalhost = [[] for e in indxsersfgrd]
        fluxhost = [[] for e in indxsersfgrd]
        if gdat.numbener > 1:
            sindhost = [[] for e in indxsersfgrd]
        sizehost = [[] for e in indxsersfgrd]
        for e in indxsersfgrd:
            lgalhost[e] = samp[getattr(gdat, strgmodl + 'indxfixplgalhostisf%d' % e)]
            bgalhost[e] = samp[getattr(gdat, strgmodl + 'indxfixpbgalhostisf%d' % e)]
            fluxhost[e] = samp[getattr(gdat, strgmodl + 'indxfixpfluxhostisf%d' % e)]
            if gdat.numbener > 1:
                sindhost[e] = samp[getattr(gdat, strgmodl + 'indxfixpsindhostisf%d' % e)]
            sizehost[e] = samp[getattr(gdat, strgmodl + 'indxfixpsizehostisf%d' % e)]
    if lensmodltype != 'none':
        beinhost = [[] for e in indxsersfgrd]
        for e in indxsersfgrd:
            if raww:
                beinhost[e] = 0.
            else:
                beinhost[e] = samp[getattr(gdat, strgmodl + 'indxfixpbeinhostisf%d' % e)]
    if hostemistype != 'none':
        ellphost = [[] for e in indxsersfgrd]
        anglhost = [[] for e in indxsersfgrd]
        serihost = [[] for e in indxsersfgrd]
        for e in indxsersfgrd:
            ellphost[e] = samp[getattr(gdat, strgmodl + 'indxfixpellphostisf%d' % e)]
            anglhost[e] = samp[getattr(gdat, strgmodl + 'indxfixpanglhostisf%d' % e)]
            serihost[e] = samp[getattr(gdat, strgmodl + 'indxfixpserihostisf%d' % e)]
    if lensmodltype != 'none':
        numbpixltemp = gdat.numbpixlcart
        defl = np.zeros((numbpixltemp, 2))
        
    # prepare the evaluation np.arrays
    if numbtrap > 0:
        initchro(gdat, gdatmodi, 'evalelem')
        if gdat.verbtype > 2:
            print 'Constructing evaluation elements...'
        
        # fill the evaluation dictionary
        indxgrideval = gdat.indxgrid
        indxpopleval = indxpopl
        indxelemeval = [[] for l in indxpopl]
        numbelem = numbelem
        # common dictionary
        for l in indxpopl:
            indxelemeval[l] = np.arange(numbelem[l])
            for strgfeat in liststrgfeateval[l]:
                if strgfeat == 'spec':
                    if boolelemlghtspat[l]:
                        sindcolr = [dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                        dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], sind=dictelem[l]['sind'], curv=dictelem[l]['curv'], \
                                  expc=dictelem[l]['expc'], sindcolr=sindcolr, spectype=spectype[l])
                    if elemtype[l].startswith('lghtline'):
                        if elemtype[l] == 'lghtlinevoig':
                            dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], elin=dictelem[l]['elin'], sigm=dictelem[l]['sigm'], \
                                                                                                            gamm=dictelem[l]['gamm'], spectype=spectype[l])
                        else:
                            dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], elin=dictelem[l]['elin'], edisintp=gdat.edisintp, spectype=spectype[l])
                else:
                    dictelem[l][strgfeat] = dictelem[l][strgfeat]

        if gdat.verbtype > 2:
            for l in indxpopl:
                print 'l'
                print l
                for strgfeat in liststrgfeateval[l]:
                    print strgfeat
                    print dictelem[l][strgfeat]
        stopchro(gdat, gdatmodi, 'evalelem')
    else:
        indxgrideval = gdat.indxgrid
        #numbelem = fittnumbelemzero
    
    # determine the indices of the pixels over which element kernels will be evaluated
    if gdat.numbpixlfull > 1:
        if numbtrap > 0:
            listindxpixlelem = [[] for l in indxpopl]
            listindxpixlelemconc = [[] for l in indxpopl]
            for l in indxpopl:
                if numbelem[l] > 0:
                    listindxpixlelem[l], listindxpixlelemconc[l] = retr_indxpixlelemconc(gdat, strgmodl, dictelem, l)
    
    if lensmodltype != 'none':
        if raww:
            sherextr = 0
        else:
            sherextr = samp[getattr(gdat, strgmodl + 'indxfixpsherextr')]
        sangextr = samp[getattr(gdat, strgmodl + 'indxfixpsangextr')]
       
        ## host halo deflection
        initchro(gdat, gdatmodi, 'deflhost')
        deflhost = [[] for e in indxsersfgrd]
            
        indxpixlmiss = gdat.indxpixlcart

        for e in indxsersfgrd:
            if gdat.verbtype > 2:
                print 'Evaluating the deflection field due to host galaxy %d' % e
                print 'lgalhost[e]'
                print lgalhost[e]
                print 'bgalhost[e]'
                print bgalhost[e]
                print 'beinhost[e]'
                print beinhost[e]
                print 'ellphost[e]'
                print ellphost[e]
                print 'anglhost[e]'
                print anglhost[e]
                print

            deflhost[e] = retr_defl(gdat, indxpixlmiss, lgalhost[e], bgalhost[e], beinhost[e], ellp=ellphost[e], angl=anglhost[e])
             
            if gdat.booldiagmode:
                indxpixltemp = slice(None)
            
            setattr(gdatobjt, strgpfix + 'deflhostisf%d' % e, deflhost[e])
       
            if gdat.verbtype > 2:
                print 'deflhost[e]'
                summgene(deflhost[e])
                
            defl += deflhost[e]
            if gdat.verbtype > 2:
                print 'After adding the host deflection...'
                print 'defl'
                summgene(defl)
        if gdat.booldiagmode:
            if not np.isfinite(deflhost).all():
                raise Exception('')
        
        stopchro(gdat, gdatmodi, 'deflhost')

        ## external shear
        initchro(gdat, gdatmodi, 'deflextr')
        deflextr = []
        indxpixltemp = gdat.indxpixlcart
        deflextr = retr_deflextr(gdat, indxpixltemp, sherextr, sangextr)
        defl += deflextr
        if gdat.verbtype > 2:
            print 'After adding the external deflection...'
            print 'defl'
            summgene(defl)
        stopchro(gdat, gdatmodi, 'deflextr')
    
    ## construct the PSF to be convolved with the image
    if gdat.pixltype == 'cart' and (psfnevaltype == 'conv' or psfnevaltype == 'full'):
        initchro(gdat, gdatmodi, 'psfnconv')
        if gdat.verbtype > 2:
            print 'Evaluating the PSF convolution kernel...'
        psfnconv = [[[] for i in gdat.indxener] for m in gdat.indxevtt]
        if gdat.pixltype == 'cart':
            if psfntype != 'singgaus':
                psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, strgmodl)
                fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            for mm, m in enumerate(gdat.indxevtt):
                for ii, i in enumerate(gdat.indxener):
                    if psfntype == 'singgaus':
                        sigm = psfp[i+m*gdat.numbener]
                    else:
                        sigm = fwhm[i, m] / 2.355
                    psfnconv[mm][ii] = AiryDisk2DKernel(sigm / gdat.sizepixl)
        
        # this is needed for state updates
        #setattr(gdatobjt, strgpfix + 'psfp', psfp)
        
        setattr(gdatobjt, strgpfix + 'psfnconv', psfnconv)
        stopchro(gdat, gdatmodi, 'psfnconv')
    
    if (psfnevaltype == 'kern' or psfnevaltype == 'full') and numbtrap > 0:
        
        psfntype = getattr(gdat, strgmodl + 'psfntype')
        if strgmodl == 'fitt':
            boolmodipsfn = getattr(gdat, strgmodl + 'boolmodipsfn')
        
        if boolinit or boolmodipsfn and strgstat == 'next':
            
            if gdat.pixltype == 'heal':
                psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, strgmodl)
                psfnintp = interp1d_pick(gdat.binsangl, psfn, axis=1)
                fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            if gdat.pixltype == 'cart':
                if gdat.kernevaltype == 'ulip':
                    psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, strgmodl)
                    psfnintp = interp1d_pick(gdat.binsangl, psfn, axis=1)

                if gdat.kernevaltype == 'bspx':
                    
                    psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsanglcart.flatten(), psfntype, strgmodl)
                    
                    # side length of the upsampled kernel
                    gdat.numbsidekernusam = 100
                    # side length of the original kernel
                    gdat.numbsidekern = gdat.numbsidekernusam / factkernusam 
                    gdat.indxsidekern = np.arange(gdat.numbsidekern)

    	        	# pad by one row and one column
    	        	#psf = np.zeros((gdat.numbsidekernusam+1, gdat.numbsidekernusam+1))
    	        	#psf[0:gdat.numbsidekernusam, 0:gdat.numbsidekernusam] = psf0
		        	
    	        	# make design matrix for each factkernusam x factkernusam region
                    nx = factkernusam + 1
                    y, x = mgrid[0:nx, 0:nx] / float(factkernusam)
                    x = x.flatten()
                    y = y.flatten()
                    print 'x'
                    print x
                    kernmatrdesi = np.array([full(nx*nx, 1), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y]).T
    	        	
                    # output np.array of coefficients
                    psfnintp = np.empty((gdat.numbsidekern, gdat.numbsidekern, kernmatrdesi.shape[1]))

    	        	# solve p = kernmatrdesi psfnintp for psfnintp
                    for iy in gdat.indxsidekern:
                        for ix in gdat.indxsidekern:
    	        	        p = psf[iy*factkernusam:(iy+1)*factkernusam+1, ix*factkernusam:(ix+1)*factkernusam+1].flatten()
    	        	        psfnintp[iy, ix, :] = dot(linalg.inv(dot(kernmatrdesi.T, kernmatrdesi)), dot(kernmatrdesi.T, p))
                setattr(gdatobjt, strgpfix + 'psfnintp', psfnintp)
        else:
            psfnintp = getattr(gdatobjt, strgpfixthis + 'psfnintp')
    
    sbrt = dict()
    for name in listnamediff:
        sbrt[name] = []
        
    if numbtrap > 0:
        if boolelemsbrtdfncanyy:
            sbrtdfnc = []
        if boolelemsbrtextsbgrdanyy:
            sbrtextsbgrd = []
        if boolelemdeflsubhanyy:
            deflsubh = []
        # retrieve or initialize state variable
        if boolelemsbrtdfncanyy:
            sbrtdfnc = np.zeros_like(gdat.expo)
        if boolelemdeflsubhanyy:
            deflsubh = np.zeros((gdat.numbpixl, 2))
        if boolelemsbrtextsbgrdanyy: 
            sbrtextsbgrd = np.zeros_like(gdat.expo)
        
        # element kernel evaluation
        if boolelemsbrtdfncanyy:
            initchro(gdat, gdatmodi, 'elemsbrtdfnc')
            sbrt['dfnc'] = []
            for l in indxpopl:
                if boolelemsbrtdfnc[l]:
                    for k in range(numbelem[l]):
                        if boolelemlght[l]:
                            varbamplextd = dictelem[l]['spec'][:, k]
                        if elemtype[l].startswith('clus'):
                            varbamplextd = dictelem[l]['nobj'][None, k]
                        if elemtype[l] == 'clusvari':
                            sbrtdfnc[0, listindxpixlelem[l][k], 0] += dictelem[l]['nobj'][k] / 2. / np.pi / dictelem[l]['gwdt'][k]**2 * \
                                np.exp(-0.5 * ((dictelem[l]['lgal'][k] - gdat.lgalgrid[listindxpixlelem[l][k]])**2 + \
                                            (dictelem[l]['bgal'][k] - gdat.bgalgrid[listindxpixlelem[l][k]])**2) / dictelem[l]['gwdt'][k]**2)
                            
                        if boolelempsfn[l]:
                            sbrtdfnc[:, listindxpixlelem[l][k], :] += retr_sbrtpnts(gdat, dictelem[l]['lgal'][k], \
                                                                           dictelem[l]['bgal'][k], varbamplextd, psfnintp, listindxpixlelem[l][k])
                        
                        if elemtype[l].startswith('lghtline'):
                            sbrtdfnc[:, 0, 0] += dictelem[l]['spec'][:, k]
                        
            sbrt['dfnc'] = sbrtdfnc
            
            setattr(gdatobjt, strgpfix + 'sbrtdfnc', sbrt['dfnc'])

            if gdat.booldiagmode:
                cntppntschec = retr_cntp(gdat, sbrt['dfnc'])
                numbelemtemp = 0
                for l in indxpopl:
                    if boolelemsbrtdfnc[l]:
                        numbelemtemp += np.sum(numbelem[l])
                if np.amin(cntppntschec) < -0.1:
                    print 'Diagnostic check failed.'
                    print 'cntppntschec'
                    summgene(cntppntschec)
                    print 'Point source spectral surface brightness is not positive-definite.'
                    #raise Exception('Point source spectral surface brightness is not positive-definite.')
            
            stopchro(gdat, gdatmodi, 'elemsbrtdfnc')
        
        if boolelemdeflsubhanyy:
            initchro(gdat, gdatmodi, 'elemdeflsubh')
            if raww:
                deflsubh = np.zeros_like(defl)
            else:
                if gdat.verbtype > 2:
                    print 'Perturbing subhalo deflection field'
                for l in indxpopl:
                    if elemtype[l] == 'lens':
                        for kk, k in enumerate(indxelemeval[l]):
                            asca = dictelem[l]['asca'][k]
                            acut = dictelem[l]['acut'][k]
                            if elemspatevaltype[l] == 'locl':
                                indxpixl = listindxpixlelem[l][kk]
                            else:
                                indxpixl = gdat.indxpixl
                            deflsubh[indxpixl, :] += retr_defl(gdat, indxpixl, \
                                                         dictelem[l]['lgal'][kk], dictelem[l]['bgal'][kk], dictelem[l]['defs'][kk], \
                                                         asca=asca, acut=acut)
                
                        # temp -- find out what is causing the features in the element convergence maps
                        #for kk, k in enumerate(indxelemeval[l]):
                        #    indxpixlpnts = retr_indxpixl(gdat, dictelem[l]['bgal'][kk], dictelem[l]['lgal'][kk])
                        #    if deflsubh[listindxpixlelem[l][kk], :]
                
                if gdat.verbtype > 2:
                    print 'deflsubh'
                    summgene(deflsubh)
                setattr(gdatobjt, strgpfix + 'deflsubh', deflsubh)
                
                if gdat.booldiagmode:
                    if not np.isfinite(deflsubh).all():
                        raise Exception('Element deflection is not finite.')

            defl += deflsubh
            if gdat.verbtype > 2:
                print 'After adding subhalo deflection to the total deflection'
                print 'defl'
                summgene(defl)

            stopchro(gdat, gdatmodi, 'elemdeflsubh')

        if boolelemsbrtextsbgrdanyy:
            initchro(gdat, gdatmodi, 'elemsbrtextsbgrd')
            if strgstat == 'this':
                for l in indxpopl:
                    if elemtype[l] == 'lghtgausbgrd':
                        for k in range(numbelem[l]):
                            sbrtextsbgrd[:, listindxpixlelem[l][k], :] += dictelem[l]['spec'][:, k, None, None] / \
                                    2. / np.pi / dictelem[l]['gwdt'][k]**2 * \
                                    np.exp(-0.5 * ((dictelem[l]['lgal'][k] - gdat.lgalgrid[None, listindxpixlelem[l][k], None])**2 + \
                                                (dictelem[l]['bgal'][k] - gdat.bgalgrid[None, listindxpixlelem[l][k], None])**2) / dictelem[l]['gwdt'][k]**2)
                
                setattr(gdatobjt, strgpfix + 'sbrtextsbgrd', sbrtextsbgrd)
            sbrt['extsbgrd'] = []
            sbrt['extsbgrd'] = sbrtextsbgrd
            
            if gdat.booldiagmode:
                cntppntschec = retr_cntp(gdat, sbrt['extsbgrd'])
                if np.amin(cntppntschec) < -0.1:
                    print 'Diagnostic check failed.'
                    print 'cntppntschec'
                    summgene(cntppntschec)
                    raise Exception('Point source spectral surface brightness is not positive-definite.')
        
            stopchro(gdat, gdatmodi, 'elemsbrtextsbgrd')
    
        if gdat.verbtype > 2:
            print 'Element related state variables after perturbations...'
            if boolelemsbrtdfncanyy:
                print 'sbrtdfnc'
                summgene(sbrtdfnc)
            if boolelemdeflsubhanyy:
                print 'deflsubh'
                summgene(deflsubh)
            if boolelemsbrtextsbgrdanyy:
                print 'sbrtextsbgrd'
                summgene(sbrtextsbgrd)
    
    if lensmodltype != 'none':
        
        # lensed surface brightness
        initchro(gdat, gdatmodi, 'sbrtlens')
        
        if gdat.verbtype > 2:
            print 'Evaluating lensed surface brightness...'
        
        if strgstat == 'this' or numbtrap > 0 and boolelemsbrtextsbgrdanyy:
            sbrt['bgrd'] = []
        if numbtrap > 0 and boolelemsbrtextsbgrdanyy:
            sbrt['bgrdgalx'] = []
        
        if gdat.numbener > 1:
            specsour = retr_spec(gdat, np.array([fluxsour]), sind=np.array([sindsour]))
            if gdat.verbtype > 2:
                print 'sindsour'
                print sindsour
        else:
            specsour = np.array([fluxsour])
        
        if gdat.verbtype > 2:
            print 'lgalsour'
            print lgalsour
            print 'bgalsour'
            print bgalsour
            print 'sizesour'
            print sizesour
            print 'ellpsour'
            print ellpsour
            print 'anglsour'
            print anglsour
            print 'fluxsour'
            print fluxsour
            print 'specsour'
            print specsour

        if numbtrap > 0 and boolelemsbrtextsbgrdanyy:
        
            if gdat.verbtype > 2:
                print 'Interpolating the background emission...'

            sbrt['bgrdgalx'] = retr_sbrtsers(gdat, gdat.lgalgrid[indxpixlelem[0]], gdat.bgalgrid[indxpixlelem[0]], \
                                                                            lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
            if gdat.verbtype > 2:
                print 'sbrt[bgrdgalx]'
                summgene(sbrt['bgrdgalx'])
                print 'sbrtextsbgrd'
                summgene(sbrtextsbgrd)
            sbrt['bgrd'] = sbrt['bgrdgalx'] + sbrtextsbgrd
        
            sbrt['lens'] = np.empty_like(gdat.cntpdata)
            for ii, i in enumerate(gdat.indxener):
                for mm, m in enumerate(gdat.indxevtt):
                    sbrtbgrdobjt = sp.interpolate.RectBivariateSpline(gdat.meanbgalcart, gdat.meanlgalcart, \
                                                            sbrt['bgrd'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)).T)
                    
                    bgalprim = gdat.bgalgrid[indxpixlelem[0]] - defl[indxpixlelem[0], 1]
                    lgalprim = gdat.lgalgrid[indxpixlelem[0]] - defl[indxpixlelem[0], 0]
                    # temp -- T?
                    sbrt['lens'][ii, :, m] = sbrtbgrdobjt(bgalprim, lgalprim, grid=False).flatten()
        else:
            if gdat.verbtype > 2:
                print 'Not interpolating the background emission...'
            
            sbrt['lens'] = retr_sbrtsers(gdat, gdat.lgalgrid - defl[gdat.indxpixl, 0], \
                                                   gdat.bgalgrid - defl[gdat.indxpixl, 1], \
                                                   lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
            sbrt['bgrd'] = retr_sbrtsers(gdat, gdat.lgalgrid, \
                                                   gdat.bgalgrid, \
                                                   lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
        setattr(gdatobjt, strgpfixthis + 'sbrtlens', sbrt['lens'])

        if gdat.booldiagmode:
            if not np.isfinite(sbrt['lens']).all():
                raise Exception('Lensed emission is not finite.')
            if (sbrt['lens'] == 0).all():
                raise Exception('Lensed emission is zero everynp.where.')

        stopchro(gdat, gdatmodi, 'sbrtlens')
        
    ### background surface brightness
    numbback = getattr(gdat, strgmodl + 'numbback')
    indxback = getattr(gdat, strgmodl + 'indxback')

    sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
    indxbacpback = getattr(gdat, strgmodl + 'indxbacpback')
    # temp
    #smthback = getattr(gdat, strgmodl + 'smthback')
    specback = getattr(gdat, strgmodl + 'specback')
    unifback = getattr(gdat, strgmodl + 'unifback')
    sbrtback = []
    # temp
    #sbrtback = np.empty((numbback, gdat.numbener, indxpixlelem[yy].size, gdat.numbevtt))
    
    # evaluate host galaxy surface brightness
    if hostemistype != 'none':
        initchro(gdat, gdatmodi, 'sbrthost')
        for e in indxsersfgrd:
            if gdat.verbtype > 2:
                print 'Evaluating the host galaxy surface brightness...'
            if gdat.numbener > 1:
                spechost = retr_spec(gdat, np.array([fluxhost[e]]), sind=np.array([sindhost[e]]))
            else:
                spechost = np.array([fluxhost[e]])
            
            if gdat.verbtype > 2:
                print 'lgalhost[e]'
                print lgalhost[e] * gdat.anglfact
                print 'bgalhost[e]'
                print bgalhost[e] * gdat.anglfact
                print 'spechost'
                print spechost
                print 'sizehost[e]'
                print sizehost[e]
                print 'ellphost[e]'
                print ellphost[e]
                print 'anglhost[e]'
                print anglhost[e]
                print 'serihost[e]'
                print serihost[e]
            
            sbrt['hostisf%d' % e] = retr_sbrtsers(gdat, gdat.lgalgrid, gdat.bgalgrid, lgalhost[e], \
                                                                        bgalhost[e], spechost, sizehost[e], ellphost[e], anglhost[e], serihost[e])
            
            setattr(gdatobjt, strgpfix + 'sbrthostisf%d' % e, sbrt['hostisf%d' % e])
                
        #sbrthost = sbrt['host']
        if gdat.verbtype > 2:
            for e in indxsersfgrd:
                print 'e'
                print e
                print 'sbrt[hostisf%d]'
                summgene(sbrt['hostisf%d' % e])
                print
        stopchro(gdat, gdatmodi, 'sbrthost')
    
    ## model emission
    initchro(gdat, gdatmodi, 'sbrtmodl')
    if gdat.verbtype > 2:
        print 'Summing up the model emission...'
    
    sbrt['modlraww'] = np.zeros((gdat.numbener, gdat.numbpixlcart, gdat.numbevtt))
    for name in listnamediff:
        if name.startswith('back'):
            indxbacktemp = int(name[4:8])
            
            if gdat.pixltype == 'heal' and (psfnevaltype == 'full' or psfnevaltype == 'conv') and not unifback[indxbacktemp]:
                sbrttemp = getattr(gdat, strgmodl + 'sbrtbackhealfull')[indxbacktemp]
            else:
                sbrttemp = sbrtbacknorm[indxbacktemp]
           
            if specback[indxbacktemp]:
                sbrt[name] = sbrttemp * bacp[indxbacpback[indxbacktemp]]
            else:
                sbrt[name] = sbrttemp * bacp[indxbacpback[indxbacktemp][gdat.indxener]][:, None, None]
        
        sbrt['modlraww'] += sbrt[name]
        if gdat.booldiagmode:
            if np.amax(sbrttemp) == 0.:
                print 'name'
                print name
                raise Exception('')

        if gdat.verbtype > 2:
            print 'name'
            print name
            print 'sbrt[name]'
            summgene(sbrt[name])
    if gdat.verbtype > 2:
        for ii, i in enumerate(gdat.indxener):
            print 'ii, i'
            print ii, i
            for mm, m in enumerate(gdat.indxevtt):
                print 'mm, m'
                print mm, m
                print 'sbrt[modlraww][ii, :, mm]'
                summgene(sbrt['modlraww'][ii, :, mm])
        print
    
    # convolve the model with the PSF
    if convdiffanyy and (psfnevaltype == 'full' or psfnevaltype == 'conv'):
        sbrt['modlconv'] = []
        # temp -- isotropic background proposals are unnecessarily entering this clause
        if gdat.verbtype > 2:
            print 'Convolving the model image with the PSF...' 
        sbrt['modlconv'] = np.zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for ii, i in enumerate(gdat.indxener):
            for mm, m in enumerate(gdat.indxevtt):
                if gdat.strgcnfg == 'pcat_ferm_igal_mock_test':
                    print 'Convolving ii, i, mm, m'
                    print ii, i, mm, m
                    print
                if gdat.pixltype == 'cart':
                    if gdat.numbpixl == gdat.numbpixlcart:
                        sbrt['modlconv'][ii, :, mm] = convolve_fft(sbrt['modlraww'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)), \
                                                                                                                             psfnconv[mm][ii]).flatten()
                    else:
                        sbrtfull = np.zeros(gdat.numbpixlcart)
                        sbrtfull[gdat.indxpixlrofi] = sbrt['modlraww'][ii, :, mm]
                        sbrtfull = sbrtfull.reshape((gdat.numbsidecart, gdat.numbsidecart))
                        sbrt['modlconv'][ii, :, mm] = convolve_fft(sbrtfull, psfnconv[mm][ii]).flatten()[gdat.indxpixlrofi]
                    indx = np.where(sbrt['modlconv'][ii, :, mm] < 1e-50)
                    sbrt['modlconv'][ii, indx, mm] = 1e-50
                if gdat.pixltype == 'heal':
                    sbrt['modlconv'][ii, :, mm] = hp.smoothing(sbrt['modlraww'][ii, :, mm], fwhm=fwhm[i, m])[gdat.indxpixlrofi]
                    sbrt['modlconv'][ii, :, mm][np.where(sbrt['modlraww'][ii, :, mm] <= 1e-50)] = 1e-50
        
        setattr(gdatobjt, strgpfix + 'sbrtmodlconv', sbrt['modlconv'])
        # temp -- this could be made faster -- need the copy() statement because sbrtdfnc gets added to sbrtmodl afterwards
        sbrt['modl'] = np.copy(sbrt['modlconv'])
    else:
        if gdat.verbtype > 2:
            print 'Skipping PSF convolution of the model...'
        sbrt['modl'] = np.copy(sbrt['modlraww'])
    
    if gdat.verbtype > 2:
        print 'sbrt[modl]'
        summgene(sbrt['modl'])
        print

    ## add PSF-convolved delta functions to the model
    if numbtrap > 0 and boolelemsbrtdfncanyy:
        if gdat.verbtype > 2:
            print 'Adding delta functions into the model...'
            print 'sbrt[dfnc]'
            summgene(sbrt['dfnc'])
            print
        sbrt['modl'] += sbrt['dfnc']
    stopchro(gdat, gdatmodi, 'sbrtmodl')
    
    if gdat.verbtype > 2:
        print 'sbrt[modl]'
        summgene(sbrt['modl'])
        print

    ### count map
    initchro(gdat, gdatmodi, 'expo')
    cntp = dict()
    cntp['modl'] = retr_cntp(gdat, sbrt['modl'])
    
    if gdat.booldiagmode:
        setattr(gdatobjt, strgpfix + 'cntpmodl', cntp['modl'])
    stopchro(gdat, gdatmodi, 'expo')

    # mock data specific
    if strgmodl == 'true' and strgstat == 'this':
        
        if raww:
            strgvarb = 'truecntpmodlraww'
        else:
            strgvarb = 'cntpdata'
        
        # generate count data
        cntptemp = []
        cntptemp = np.zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for i in gdat.indxener:
            for j in gdat.indxpixl:
                for m in gdat.indxevtt:
                    cntptemp[i, j, m] = poisson(cntp['modl'][i, j, m])
        setattr(gdat, strgvarb + '', cntptemp)
        setattr(gdat, strgvarb, cntptemp)
    
        if not gdat.boolsqzeexpo and np.amax(cntptemp) == 0 and not raww:
            print 'gdat.deltener'
            summgene(gdat.deltener)
            print 'gdat.apix'
            print gdat.apix
            print gdat.apix * gdat.anglfact**2
            print 'gdat.expo'
            summgene(gdat.expo)
            print 'gdat.cntpdata'
            summgene(gdat.cntpdata)
            print 'cntp[modl]'
            summgene(cntp['modl'])
            raise Exception('Data is zero.')
        
        if raww:
            return
        else:
            retr_datatick(gdat)
    
        proc_cntpdata(gdat)
    
    ## diagnostics
    if gdat.booldiagmode:
        frac = cntp['modl'] / np.mean(cntp['modl'])
        if np.amin(frac) < -1e-3 and np.amin(cntp['modl']) < -0.1:
            print 'Total model surface brightness is not positive-definite.'
            print 'cntp[modl]'
            summgene(cntp['modl'])
            print 'frac'
            summgene(frac)
            raise Exception('')
        
        indxcubebadd = np.where(cntp['modl'] < 0.)[0]
        if indxcubebadd.size > 0:
            print 'Warning! Model prediction is negative. Correcting to 1e-20...'
            cntp['modl'][indxcubebadd] = 1e-20
    stopchro(gdat, gdatmodi, 'modl')

    # log-prior
    initchro(gdat, gdatmodi, 'lpri')
    if gdat.verbtype > 2:
        print 'Evaluating the prior...'
        
    lpri = np.zeros(numblpri)
    if numbtrap > 0:
        
        for l in indxpopl:
            lpri[0] -= 0.5 * gdat.priofactdoff * numbcomp[l] * numbelem[l]
        
        if gdat.penalpridiff:
            sbrtdatapnts = gdat.sbrtdata - sbrt['dfnc']
            if gdat.pixltype == 'heal':
                raise Exception('')
            if gdat.pixltype == 'cart':
                psecodimdatapnts = np.empty((gdat.numbener, gdat.numbsidecart / 2, gdat.numbevtt))
                psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, strgmodl)
                fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
                sigm = fwhm / 2.355
                psecodimdatapntsprio = np.exp(-2. * gdat.meanmpolodim[None, :, None] / (0.1 / sigm[:, None, :]))
                lpridiff = 0.
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        psecdatapnts = retr_psec(gdat, sbrtdatapnts[i, :, m])
                        psecodimdatapnts[i, :, m] = retr_psecodim(gdat, psecdatapnts)
                        psecodimdatapnts[i, :, m] /= psecodimdatapnts[i, 0, m]
                        lpridiff += -0.5 * np.sum((psecodimdatapnts[i, :, m] - psecodimdatapntsprio[i, :, m])**2)
                        setattr(gdatobjt, strgpfix + 'psecodimdatapntsen%02devt%d' % (i, m), psecodimdatapnts[i, :, m])
                        setattr(gdatobjt, strgpfix + 'psecodimdatapntsprioen%02devt%d'% (i, m), psecodimdatapntsprio[i, :, m])
            lpri[1] = lpridiff 
            setattr(gdatobjt, strgpfix + 'lpridiff', lpridiff)
        
        meanelem = samp[indxfixpmeanelem]
        for l in indxpopl:
            lpri[2] += retr_lprbpois(numbelem[l], meanelem[l])
        
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
        for l in indxpopl:
            for g, (strgfeat, strgpdfn) in enumerate(zip(liststrgfeatprio[l], liststrgpdfnprio[l])):
                indxlpritemp = 3 + l * maxmnumbcomp + g
                lpri[indxlpritemp] = retr_lprielem(gdat, strgmodl, l, g, strgfeat, strgpdfn, samp, dictelem, numbelem)
    lpritotl = np.sum(lpri)
    
    if gdat.verbtype > 1:
        print 'lpritotl'
        print lpritotl
    
    ### log-likelihood
    initchro(gdat, gdatmodi, 'llik')
    llik = retr_llik(gdat, strgmodl, cntp['modl'])
    
    if gdat.verbtype > 2:
        print 'cntp[modl]'
        summgene(cntp['modl'])
        print 'np.sum(cntp[modl], (1, 2))'
        print np.sum(cntp['modl'], (1, 2))
        print 'np.sum(gdat.cntpdata, (1, 2))'
        print np.sum(gdat.cntpdata, (1, 2))

    if gdat.booldiagmode:
        if not np.isfinite(llik).all():
            raise Exception('Likelihood is not finite.')
    
    lliktotl = np.sum(llik)
    if gdat.booldiagmode:
        if isinstance(lliktotl, np.ndarray):
            raise Exception('')
        if not np.isfinite(lliktotl).all():
            raise Exception('')

    numbfixp = getattr(gdat, strgmodl + 'numbfixp')
    numbdoff = gdat.numbdata - numbfixp
    if numbtrap > 0:
        for l in indxpopl:
            numbdoff -= len(indxsampcomp['comp'][l])

    setattr(gdatobjt, strgpfix + 'llik', llik)
    setattr(gdatobjt, strgpfix + 'llikmean', lliktotl / gdat.numbdata) 
    setattr(gdatobjt, strgpfix + 'llikcmea', lliktotl / (gdat.numbdata - numbdoff)) 

    if gdat.verbtype > 2:
        print 'llik'
        summgene(llik)
    if gdat.verbtype > 1:
        print 'lliktotl'
        print lliktotl
    stopchro(gdat, gdatmodi, 'llik')

    lpostotl = lpritotl + lliktotl
    if gdat.verbtype > 1:
        print 'lpostotl'
        print lpostotl

    setattr(gdatobjt, strgpfix + 'lpritotl', lpritotl) 
    setattr(gdatobjt, strgpfix + 'lliktotl', lliktotl)
    setattr(gdatobjt, strgpfix + 'lpostotl', lpostotl) 
    
    stopchro(gdat, gdatmodi, 'lpri')
    
    if strgstat == 'next':
        return

    initchro(gdat, gdatmodi, 'tert')
    
    setattr(gdatobjt, strgpfix + 'lpri', lpri)
    
    if numbtrap > 0:
        setattr(gdatobjt, strgpfix + 'lpripena', lpri[0])
    
    dicttert = {}
    
    if strgstat == 'this' and numbtrap > 0:
        numbelempopl = np.zeros(numbpopl)
        for l in indxpopl:
            numbelempopl[l] = np.sum(numbelem[l])
        setattr(gdatobjt, strgpfix + 'numbelem', numbelem)
        setattr(gdatobjt, strgpfix + 'numbelempopl', numbelempopl)
    
    ## load necessary variables
        
    numbdeflsubhplot = getattr(gdat, strgmodl + 'numbdeflsubhplot')
    numbdeflsingplot = getattr(gdat, strgmodl + 'numbdeflsingplot')
    liststrgfeatcorr = getattr(gdat, strgmodl + 'liststrgfeatcorr')
    liststrgfeatcorrtotl = getattr(gdat, strgmodl + 'liststrgfeatcorrtotl')
    numbback = getattr(gdat, strgmodl + 'numbback')
    nameback = getattr(gdat, strgmodl + 'nameback')
    numblablsbrt = getattr(gdat, strgmodl + 'numblablsbrt')
    listnameecomtotl = getattr(gdat, strgmodl + 'listnameecomtotl')
    listnamegcom = getattr(gdat, strgmodl + 'listnamegcom')
    
    ## set previously calculated tertiary variables
    setattr(gdatobjt, strgpfix + 'bacp', bacp)
    if psfnevaltype != 'none':
        setattr(gdatobjt, strgpfix + 'psfp', psfp)
    
    ## derived variables
    ## residual count map 
    cntp['resi'] = []
    cntp['resi'] = gdat.cntpdata - cntp['modl']
    
    setattr(gdatobjt, strgpfix + 'cntpmodl', cntp['modl'])
    setattr(gdatobjt, strgpfix + 'cntpresi', cntp['resi'])
    setattr(gdatobjt, strgpfix + 'llik', llik)
    #if lensmodltype != 'none':
    #    setattr(gdatobjt, strgpfix + 'deflhost', deflhost)
    
    namefeatsort = getattr(gdat, strgmodl + 'namefeatsort')
                      
    if lensmodltype != 'none' or hostemistype != 'none':
        numbsersfgrd = getattr(gdat, strgmodl + 'numbsersfgrd')
    if lensmodltype != 'none':
        
        setattr(gdatobjt, strgpfix + 'defl', defl)
        massfrombein = getattr(gdat, strgmodl + 'massfrombein')
        for e in indxsersfgrd:
            masshostbein = massfrombein * beinhost[e]**2
            setattr(gdatobjt, strgpfix + 'masshostisf%dbein' % e, masshostbein)
        ### sort with respect to deflection at scale radius
        if numbtrap > 0:
            for l in indxpopl:
                if numbelem[l] > 0:
                    indxelemsortampl = np.argsort(dictelem[l][namefeatsort[l]])[::-1]
                    for strgcomp in liststrgcomp[l]:
                        dictelem[l][strgcomp + 'sort'] = dictelem[l][strgcomp][indxelemsortampl]

        deflsing = np.zeros((gdat.numbpixlcart, 2, numbdeflsingplot))
        conv = np.zeros((gdat.numbpixlcart))
        convpsec = np.zeros(((gdat.numbsidecart / 2)**2))
        convpsecodim = np.zeros((gdat.numbsidecart / 2))
        if numbtrap > 0:
            if boolelemlens:
                indxpopllens = elemtype.index('lens')
        numbdeflsing = 2
        if numbtrap > 0:
            if boolelemlens:
                if numbelem[indxpopllens] > 0:
                    numbdeflsing += min(numbdeflsubhplot, numbelem[indxpopllens]) 
                    numbdeflsing += 1
                for k in range(numbdeflsing):
                    indxpixltemp = gdat.indxpixlcart
                    if k == 0:
                        # temp -- should take other sersics into account
                        deflsing[indxpixltemp, :, k] = deflhost[0]
                    elif k == 1:
                        deflsing[indxpixltemp, :, k] = deflextr
                    elif k == 2:
                        deflsing[indxpixltemp, :, k] = defl - deflextr - deflhost[0]
                    else:
                        asca = dictelem[indxpopllens]['ascasort'][None, k-3]
                        acut = dictelem[indxpopllens]['acutsort'][None, k-3]
                        deflsing[listindxpixlelem[indxpopllens][k], :, k] = retr_defl(gdat, listindxpixlelem[indxpopllens][k], \
                                                      dictelem[indxpopllens]['lgalsort'][None, k-3], dictelem[indxpopllens]['bgalsort'][None, k-3], \
                                                      dictelem[indxpopllens]['defssort'][None, k-3], asca=asca, acut=acut)

        # convergence
        ## total
        conv[:] = retr_conv(gdat, defl) 
        convhost = np.zeros((numbsersfgrd, gdat.numbpixlcart))
        for e in indxsersfgrd:
            convhost[e, :] = retr_conv(gdat, deflhost[e]) 
        
        ### power spectrum
        #### two dimensional
        convpsec[:] = retr_psec(gdat, conv[:])
        
        #### one dimensional
        convpsecodim[:] = retr_psecodim(gdat, convpsec[:]) 
        setattr(gdatobjt, strgpfix + 'convpsec', convpsec)
        setattr(gdatobjt, strgpfix + 'convpsecodim', convpsecodim)
        setattr(gdatobjt, strgpfix + 'conv', conv[...])
        for e in indxsersfgrd:
            setattr(gdatobjt, strgpfix + 'convisf%d' % e, convhost[e, ...])
        
        ## subhalos
        if numbtrap > 0:
            if boolelemlens:
                convelem = np.zeros((gdat.numbpixl))
                convpsecelem = np.zeros(((gdat.numbsidecart / 2)**2))
                convpsecelemodim = np.zeros((gdat.numbsidecart / 2))
                ### convergence
                convelem[:] = retr_conv(gdat, deflsubh) 
                ###  power spectrum
                ##### two dimensional
                convpsecelem[:] = retr_psec(gdat, convelem[:])
                ##### one dimensional
                convpsecelemodim[:] = retr_psecodim(gdat, convpsecelem[:]) 
                setattr(gdatobjt, strgpfix + 'convpsecelem', convpsecelem)
                setattr(gdatobjt, strgpfix + 'convpsecelemodim', convpsecelemodim)
                setattr(gdatobjt, strgpfix + 'convelem', convelem[...])
                setattr(gdatobjt, strgpfix + 'defl', defl)
        
        ### magnification
        magn = np.empty((gdat.numbpixlcart))
        histdefl = np.empty((gdat.numbdefl))
        if numbtrap > 0 and boolelemlens:
            histdeflsubh = np.empty((gdat.numbdefl))
        deflsingmgtd = np.zeros((gdat.numbpixlcart, numbdeflsingplot))
        magn[:] = 1. / retr_invm(gdat, defl) 
        histdefl[:] = np.histogram(defl, bins=gdat.binsdefl)[0]
        if numbtrap > 0:
            if boolelemlens:
                histdeflsubh[:] = np.histogram(deflsubh, bins=gdat.binsdeflsubh)[0]
        deflsingmgtd[:, :] = np.sqrt(np.sum(deflsing[...]**2, axis=1))
        if numbtrap > 0:
            if boolelemlens:
                setattr(gdatobjt, strgpfix + 'histdeflsubh', histdeflsubh)
        setattr(gdatobjt, strgpfix + 'histdefl', histdefl)
        setattr(gdatobjt, strgpfix + 'magn', magn[...])
        setattr(gdatobjt, strgpfix + 'deflsing', deflsing[...])
        setattr(gdatobjt, strgpfix + 'deflsingmgtd', deflsingmgtd[...])
    
    liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
    liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
    liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
    liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
    
    ## element related
    if numbtrap > 0:
        if gdat.numbpixl == 1:
            for l in indxpopl:
                for k in range(numbelem[l]):
                    setattr(gdatobjt, strgpfix + 'speclinepop%d%04d' % (l, k), dictelem[l]['spec'][:, k])
        
        if gdat.datatype == 'mock' and strgmodl == 'true' and gdat.numbpixl > 1:
            gdat.refrlgal = [[] for l in gdat.trueindxpopl]
            gdat.refrbgal = [[] for l in gdat.trueindxpopl]
            for l in gdat.trueindxpopl:
                gdat.refrlgal[l] = np.tile(dictelem[l]['lgal'], [3] + list(np.ones(dictelem[l]['lgal'].ndim, dtype=int)))
                gdat.refrbgal[l] = np.tile(dictelem[l]['bgal'], [3] + list(np.ones(dictelem[l]['bgal'].ndim, dtype=int)))
    
        for l in indxpopl:
            if elemtype[l] == 'lghtpntspuls':
                dictelem[l]['per1'] = retr_per1(dictelem[l]['per0'], dictelem[l]['magf'])
        
    if numbtrap > 0:
        if strgstat == 'this' or gdat.boolrefeforc and strgmodl == 'fitt':
            # correlate the fitting model elements with the reference elements
            if gdat.refrinfo and not (strgmodl == 'true' and gdat.datatype == 'mock') and gdat.boolasscrefr:
                indxelemrefrasschits = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                indxelemfittasschits = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        if gdat.refrnumbelem[q] == 0:
                            continue
                        
                        indxelemfittmatr = np.empty((gdat.refrnumbelem[q], numbelem[l]), dtype=int)
                        indxelemrefrmatr = np.empty((gdat.refrnumbelem[q], numbelem[l]), dtype=int)
                        matrdist = np.empty((gdat.refrnumbelem[q], numbelem[l]))
                        for k in range(numbelem[l]):
                            # construct a matrix of angular distances between reference and fitting elements
                            if elemtype[l].startswith('lghtline'):
                                matrdist[:, k] = abs(gdat.refrelin[q][0, :] - dictelem[l]['elin'][k]) / gdat.refrelin[q][0, :]
                            else:
                                matrdist[:, k] = retr_angldist(gdat, gdat.refrlgal[q][0, :], gdat.refrbgal[q][0, :], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k])
                            indxelemrefrmatr[:, k] = np.arange(gdat.refrnumbelem[q])
                            indxelemfittmatr[:, k] = k
                        matrdist = matrdist.flatten()
                        indxelemrefrmatr = indxelemrefrmatr.flatten()
                        indxelemfittmatr = indxelemfittmatr.flatten()

                        # take only angular separations smaller than some threshold
                        indxmatrthrs = np.where(matrdist < gdat.anglassc)
                        matrdist = matrdist[indxmatrthrs]
                        indxelemrefrmatr = indxelemrefrmatr[indxmatrthrs]
                        indxelemfittmatr = indxelemfittmatr[indxmatrthrs]

                        # sort the remaining associations with respect to distance
                        indxmatrsort = np.argsort(matrdist)
                        matrdist = matrdist[indxmatrsort]
                        indxelemrefrmatr = indxelemrefrmatr[indxmatrsort]
                        indxelemfittmatr = indxelemfittmatr[indxmatrsort]
                        
                        for c in range(matrdist.size):
                            if indxelemrefrmatr[c] in indxelemrefrasschits[q][l] or indxelemfittmatr[c] in indxelemfittasschits[q][l]:
                                continue
                            indxelemrefrasschits[q][l].append(indxelemrefrmatr[c])
                            indxelemfittasschits[q][l].append(indxelemfittmatr[c])
                        
                        indxelemrefrasschits[q][l] = np.array(indxelemrefrasschits[q][l])
                        indxelemfittasschits[q][l] = np.array(indxelemfittasschits[q][l])
                setattr(gdatobjt, strgpfix + 'indxelemrefrasschits', indxelemrefrasschits)
                setattr(gdatobjt, strgpfix + 'indxelemfittasschits', indxelemfittasschits)
                
                indxelemrefrasscmiss = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                indxelemfittasscfals = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        # indices of the reference elements not associated with the fitting model elements
                        if gdat.refrnumbelem[q] > 0:
                            indxelemrefrasscmiss[q][l] = np.setdiff1d(np.arange(gdat.refrnumbelem[q]), indxelemrefrasschits[q][l])
                        # indices of the fitting model elements not associated with the reference elements
                        if numbelem[l] > 0:
                            indxelemfittasscfals[q][l] = np.setdiff1d(np.arange(numbelem[l]), indxelemfittasschits[q][l])
                setattr(gdatobjt, strgpfix + 'indxelemrefrasscmiss', indxelemrefrasscmiss)
                setattr(gdatobjt, strgpfix + 'indxelemfittasscfals', indxelemfittasscfals)
                
                for q in gdat.indxrefr:
                    if gdat.refrnumbelempopl[q] == 0:
                        continue
                    for l in gdat.fittindxpopl:
                        # collect the associated reference element feature for each fitting element 
                        for strgfeat in gdat.refrliststrgfeatonly[q][l]:
                            name = strgfeat + gdat.listnamerefr[q]
                            refrfeat = getattr(gdat, 'refr' + strgfeat)
                            dictelem[l][name] = np.zeros(numbelem[l])
                            if len(refrfeat[q]) > 0 and len(indxelemrefrasschits[q][l]) > 0:
                                dictelem[l][name][indxelemfittasschits[q][l]] = refrfeat[q][0, indxelemrefrasschits[q][l]]
                        
                        # collect the error in the associated reference element amplitude
                        for strgfeat in gdat.liststrgfeatcomm[q][l]:
                            refrfeat = getattr(gdat, 'refr' + strgfeat)
                            if strgfeat == gdat.fittnamefeatampl[l] and len(indxelemfittasschits[q][l]) > 0:
                                dictelem[l]['aerr' + gdat.listnamerefr[q]] = np.zeros(numbelem[l])
                                fittfeattemp = dictelem[l][strgfeat][indxelemfittasschits[q][l]]
                                refrfeattemp = refrfeat[q][0, indxelemrefrasschits[q][l]]
                                if gdat.booldiagmode:
                                    if not np.isfinite(refrfeattemp).all():
                                        print 'refrfeattemp'
                                        print refrfeattemp
                                        raise Exception('')
                                dictelem[l]['aerr' + gdat.listnamerefr[q]][indxelemfittasschits[q][l]] = 100. * (fittfeattemp - refrfeattemp) / refrfeattemp
                
            if gdat.boolrefeforc and strgmodl == 'fitt':
                for l in gdat.fittindxpopl:
                    for strgfeat in gdat.fittliststrgfeat[l]:
                        if strgfeat in gdat.refrliststrgfeat[gdat.indxrefrforc[l]]:
                            if len(indxelemrefrasschits[gdat.indxrefrforc[l]][l]) == 0:
                                continue
                            refrfeat = getattr(gdat, 'refr' + strgfeat)[gdat.indxrefrforc[l]][0, indxelemrefrasschits[gdat.indxrefrforc[l]][l]]
                            if len(dictelem[l][strgfeat]) == 0:
                                continue
                            lpritotl += -2. * np.sum(1e6 * (dictelem[l][strgfeat][indxelemfittasschits[gdat.indxrefrforc[l]][l]] - refrfeat)**2 / refrfeat**2)

    # other tertiary variables continues
    ## number of degrees of freedom

    chi2doff = np.sum(cntp['resi']**2 / gdat.varidata) / numbdoff
    if gdat.booldiagmode:
        if not np.isfinite(cntp['resi']).all():
            raise Exception('')
        if not np.isfinite(numbdoff):
            raise Exception('')
        if not np.isfinite(chi2doff):
            print 'numbdoff'
            print numbdoff
            print 'gdat.cntpdata'
            summgene(gdat.cntpdata)
            print 'cntp[resi]'
            summgene(cntp['resi'])
            raise Exception('')
    setattr(gdatobjt, strgpfix + 'numbdoff', numbdoff)
    setattr(gdatobjt, strgpfix + 'chi2doff', chi2doff)
    
    if numbtrap > 0:
        
        ### derived quantities
        for l in indxpopl:

            # luminosity
            if boolelemlght[l] and 'flux' in liststrgfeat[l]:
                for strgfeat in liststrgfeat[l]:
                    if strgfeat.startswith('reds') and strgfeat != 'reds':
                        namerefr = strgfeat[-4:]
                        dictelem[l]['lumi' + namerefr] = np.zeros(numbelem[l]) + np.nan
                        dictelem[l]['dlos' + namerefr] = np.zeros(numbelem[l]) + np.nan
                        reds = dictelem[l]['reds' + namerefr]
                        indxgood = np.where(np.isfinite(dictelem[l]['reds' + namerefr]))[0]
                        if indxgood.size > 0:
                            # temp -- these units only work for energy units of keV
                            dlos = gdat.adisobjt(reds)
                            dictelem[l]['dlos' + namerefr][indxgood] = dlos
                            lumi = retr_lumi(gdat, dictelem[l]['flux'], dlos, reds)
                            dictelem[l]['lumi' + namerefr][indxgood] = lumi
        
            if elemtype[l] == 'lghtpntsagnntrue':
                dictelem[l]['reds'] = gdat.redsfromdlosobjt(dictelem[l]['dlos'])
            if elemtype[l] == 'lghtpntspuls':
                dictelem[l]['mass'] = full([numbelem[l]], 3.)

            if gdat.verbtype > 2:
                print 'l'
                print l
            if gdat.numbpixl > 1:
                #### radial and angular coordinates
                dictelem[l]['gang'] = retr_gang(dictelem[l]['lgal'], dictelem[l]['bgal'])
                dictelem[l]['aang'] = retr_aang(dictelem[l]['lgal'], dictelem[l]['bgal'])
            
            if boolelemlght[l]:
                #### number of expected counts
                if gdat.numbpixlfull > 1:
                    dictelem[l]['cnts'] = retr_cntspnts(gdat, [dictelem[l]['lgal'], dictelem[l]['bgal']], dictelem[l]['spec'])
                else:
                    dictelem[l]['cnts'] = retr_cntspnts(gdat, [dictelem[l]['elin']], dictelem[l]['spec'])
            
            #### delta log-likelihood
            dictelem[l]['deltllik'] = np.zeros(numbelem[l])
            if not (strgmodl == 'true' and gdat.checprio): 
                if gdat.verbtype > 2:
                    print
                    print 'Calculating log-likelihood differences when removing elements from the model.'
                for k in range(numbelem[l]):
                    gdatmoditemp = tdpy.util.gdatstrt()
                    prep_gdatmodi(gdat, gdatmoditemp, gdatobjt, strgmodl)
                    prop_stat(gdat, gdatmoditemp, strgmodl, deth=True, thisindxpopl=l, thisindxelem=k)
                    proc_samp(gdat, gdatmoditemp, 'next', strgmodl, boolinit=True)
                    
                    if gdat.booldiagmode:
                        if not np.isfinite(lliktotl):
                            print 'lliktotl'
                            print lliktotl
                            raise Exception('')
                    
                    if strgmodl == 'true':
                        nextlliktotl = gdat.nexttruelliktotl
                    else:
                        nextlliktotl = gdatmoditemp.nextlliktotl
                    dictelem[l]['deltllik'][k] = lliktotl - nextlliktotl
                    
                    if gdat.verbtype > 2:
                        print 'k'
                        print k
                        print 'dictelem[l][deltllik][k]'
                        print dictelem[l]['deltllik'][k]
                        print
                if gdat.verbtype > 2:
                    print
                if gdat.verbtype > 2:
                    print 'deltllik calculation ended.'
                    print
    
    calcerrr = getattr(gdat, strgmodl + 'calcerrr')
    
    # more derived parameters
    if psfnevaltype == 'kern' or psfnevaltype == 'full':
        if boolinit or boolmodipsfn and strgstat == 'next':
            setattr(gdatobjt, strgpfix + 'psfn', psfn)

            ### PSF FWHM
            if gdat.pixltype == 'cart':
                fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            setattr(gdatobjt, strgpfix + 'fwhm', fwhm)
    
    if numbtrap > 0 and boolelemsbrtdfncanyy:
        boolelemsbrt = getattr(gdat, strgmodl + 'boolelemsbrt')
        
        if gdat.fittnumbtrap > 0:
            sbrt['dfnctotl'] = np.zeros_like(gdat.expo)
            sbrt['dfncsubt'] = np.zeros_like(gdat.expo)
            sbrt['dfncsupt'] = np.zeros_like(gdat.expo)
            for l in indxpopl:
                if calcerrr[l]:
                    sbrt['dfncfull'] = np.zeros_like(gdat.expo)
                if boolelemsbrt[l]:
                    for k in range(numbelem[l]):
                        
                        # read normalization from the element dictionary
                        if boolelemlght[l]:
                            varbamplextd = dictelem[l]['spec'][:, k]
                        if elemtype[l].startswith('clus'):
                            varbamplextd = dictelem[l]['nobj'][None, k]
                        
                        # calculate imprint on the element surface brightness state variable
                        if boolelempsfn[l]:
                            sbrttemp = retr_sbrtpnts(gdat, dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], \
                                                                                                    varbamplextd, psfnintp, listindxpixlelem[l][k])
                        indxpixltemp = listindxpixlelem[l][k]

                        if elemtype[l].startswith('lghtline'):
                            sbrttemp = dictelem[l]['spec'][:, k, None, None]
                        
                        # add it to the state variable depending on the significance
                        sbrt['dfnctotl'][:, indxpixltemp, :] += sbrttemp
                        if dictelem[l]['deltllik'][k] > 35:
                            sbrt['dfncsupt'][:, indxpixltemp, :] += sbrttemp
                        if dictelem[l]['deltllik'][k] < 35:
                            sbrt['dfncsubt'][:, indxpixltemp, :] += sbrttemp
                        
                        # calculate imprint without PSF truncation to calculate approximation errors
                        if calcerrr[l]:
                            sbrt['dfncfull'][:, :, :] += retr_sbrtpnts(gdat, dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], \
                                                                                                            varbamplextd, psfnintp, gdat.indxpixl)
            
                setattr(gdatobjt, strgpfix + 'sbrtdfncsubtpop%d' % l, sbrt['dfncsubt'])
                
    if numbtrap > 0 and boolelemsbrtextsbgrdanyy:
        if gdat.booldiagmode:
            numbtemp = 0
            for l in indxpopl:
                if boolelemsbrtextsbgrd[l]:
                    numbtemp += np.sum(numbelem[l])
            if numbtemp > 0 and (sbrtextsbgrd == 0.).all():
                raise Exception('')

        sbrt['bgrdexts'] = sbrtextsbgrd

    #### count maps
    if gdat.correxpo:
        cntp = dict()
        for name in listnamegcom:
            cntp[name] = retr_cntp(gdat, sbrt[name])
            setattr(gdatobjt, strgpfix + 'cntp' + name + '', cntp[name])
    
    ### spatial averages
    sbrtmean = dict()
    sbrtstdv = dict()
    for name in listnamegcom:
        sbrtmean[name], sbrtstdv[name] = retr_spatmean(gdat, sbrt[name])
        for b in gdat.indxspatmean:
            setattr(gdatobjt, strgpfix + 'sbrt%smea%d' % (name, b), sbrtmean[name][b])
            setattr(gdatobjt, strgpfix + 'sbrt%sstd%d' % (name, b), sbrtstdv[name][b])
    
    if numbtrap > 0:
        if boolelemsbrtdfncanyy:
            for i in gdat.indxener:
                if 'dark' in listnamegcom:
                    fracsdenmeandarkdfncsubt = sbrtmean['dfncsubt'][0][0][i] / (sbrtmean['dfncsubt'][0][0][i] + sbrtmean['dark'][0][0][i])
                else:
                    fracsdenmeandarkdfncsubt = 1.
                setattr(gdatobjt, strgpfix + 'fracsdenmeandarkdfncsubten%02d' % i, np.array([fracsdenmeandarkdfncsubt]))
            
            if 'dark' in listnamegcom:
                booldfncsubt = float(np.where(sbrtmean['dfncsubt'][0][0] > sbrtmean['dark'][0][0])[0].any())
            else:
                booldfncsubt = 1.
            setattr(gdatobjt, strgpfix + 'booldfncsubt', np.array([booldfncsubt]))

    # find the 1-point function of the count maps of all emission components including the total emission
    for name in listnamegcom:
        for m in gdat.indxevtt:
            if gdat.numbpixl > 1:
                for i in gdat.indxener: 
                    histcntp = np.histogram(cntp[name][i, :, m], bins=gdat.binscntpmodl)[0]
                    strgtemp = strgpfix + 'histcntp' + name + 'en%02devt%d' % (i, m)
                    setattr(gdatobjt, strgtemp, histcntp)
                    
                    if i == 0 and m == 0 and (name == 'dfnc' or name == 'dfncsubt'):
                        for strgbins in ['lowr', 'higr']:
                            strgtemp = strgpfix + 'histcntp' + strgbins + name + 'en%02devt%d' % (i, m)
                            if strgbins == 'lowr':
                                setattr(gdatobjt, strgtemp, np.array([float(np.sum(histcntp[:gdat.numbtickcbar-1]))]))
                            else:
                                setattr(gdatobjt, strgtemp, np.array([float(np.sum(histcntp[gdat.numbtickcbar-1:]))]))
            else:
                histcntp = np.histogram(cntp[name][:, 0, m], bins=gdat.binscntpmodl)[0]
                setattr(gdatobjt, strgpfix + 'histcntp' + name + 'evt%d' % m, histcntp)

    if lensmodltype != 'none':
        adishost = getattr(gdat, strgmodl + 'adishost')
        adishostsour = getattr(gdat, strgmodl + 'adishostsour')
        adissour = getattr(gdat, strgmodl + 'adissour')
        mdencrit = getattr(gdat, strgmodl + 'mdencrit')
        if strgmodl == 'true':
            s2nr = []
            s2nr = cntp['lens'] / np.sqrt(cntp['modl'])
            setattr(gdatobjt, strgpfix + 's2nr', s2nr)
        cntplensgrad = np.empty((gdat.numbener, gdat.numbpixlcart, gdat.numbevtt, 2))
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                cntplenstemp = np.zeros(gdat.numbpixlcart)
                cntplenstemp[gdat.indxpixlrofi] = cntp['lens'][i, :, m]
                cntplensgrad[i, :, m, :] = retr_gradmaps(gdat, cntplenstemp) * gdat.sizepixl
        
        cntplensgradmgtd = np.sqrt(np.sum(cntplensgrad**2, axis=3))
        cntplensgrad *= gdat.sizepixl
        indx = np.where(np.fabs(cntplensgrad) > 1. * gdat.sizepixl)
        cntplensgrad[indx] = np.sign(cntplensgrad[indx]) * 1. * gdat.sizepixl
        deflmgtd = np.sqrt(np.sum(defl**2, axis=1))
        setattr(gdatobjt, strgpfix + 'deflmgtd', deflmgtd)
        setattr(gdatobjt, strgpfix + 'cntplensgrad', cntplensgrad)
        setattr(gdatobjt, strgpfix + 'cntplensgradmgtd', cntplensgradmgtd)

    if numbtrap > 0:
        for l in indxpopl:
            if boolelemlght[l]:
                #### spectra
                if gdat.numbpixlfull > 1:
                    sindcolr = [dictelem[l]['sindcolr%04d' % i] for i in gdat.indxenerinde]
                    dictelem[l]['specplot'] = retr_spec(gdat, dictelem[l]['flux'], sind=dictelem[l]['sind'], \
                                                                 curv=dictelem[l]['curv'], expc=dictelem[l]['expc'], \
                                                                 sindcolr=sindcolr, spectype=spectype[l], plot=True)
                
                if gdat.datatype == 'inpt':
                    if gdat.exprtype == 'ferm':
                        # temp
                        try:
                            dictelem[l]['sbrt0018'] = gdat.sbrt0018objt(dictelem[l]['bgal'], dictelem[l]['lgal'])
                        except:
                            dictelem[l]['sbrt0018'] = dictelem[l]['bgal'] * 0.

            if elemtype[l] == 'lens':
                #### distance to the source
                if lensmodltype != 'none':
                    dictelem[l]['diss'] = retr_angldist(gdat, dictelem[l]['lgal'],  dictelem[l]['bgal'], lgalsour, bgalsour)
                
                if lensmodltype == 'elem' or lensmodltype == 'full':
                    dictelem[l]['deflprof'] = np.empty((gdat.numbanglfull, numbelem[l]))
                    dictelem[l]['mcut'] = np.empty(numbelem[l])
                    dictelem[l]['rele'] = np.empty(numbelem[l])
                    dictelem[l]['reln'] = np.empty(numbelem[l])
                    dictelem[l]['relk'] = np.empty(numbelem[l])
                    dictelem[l]['relf'] = np.empty(numbelem[l])
                    dictelem[l]['reld'] = np.empty(numbelem[l])
                    dictelem[l]['relc'] = np.empty(numbelem[l])
                    dictelem[l]['relm'] = np.empty(numbelem[l])

                    # temp -- this can be placed earlier in the code
                    cntplensobjt = sp.interpolate.RectBivariateSpline(gdat.meanbgalcart, gdat.meanlgalcart, \
                                                            cntp['lens'][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)).T)
                    
                    for k in np.arange(numbelem[l]):
                        
                        asca = dictelem[l]['asca'][k]
                        acut = dictelem[l]['acut'][k]
                        
                        #### deflection profiles
                        dictelem[l]['deflprof'][:, k] = retr_deflcutf(gdat.meananglfull, dictelem[l]['defs'][k], asca, acut)
         
                        ### truncated mass 
                        dictelem[l]['mcut'][k] = retr_mcut(gdat, dictelem[l]['defs'][k], asca, acut, adishost, mdencrit)

                        #### dot product with the source flux gradient
                        # temp -- weigh the energy and PSF bins
                        dictelem[l]['rele'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], \
                                                                              dictelem[l]['defs'][k], asca, acut, gdat.indxpixl)
                        
                        dictelem[l]['relf'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], \
                                                                         dictelem[l]['defs'][k], asca, acut, gdat.indxpixl, cntpmodl=cntp['modl'][0, :, 0])
                        
                        deflelem = retr_defl(gdat, gdat.indxpixl, dictelem[l]['lgal'][k], \
                                                                            dictelem[l]['bgal'][k], dictelem[l]['defs'][k], asca=asca, acut=acut)
                        bgalprim = gdat.bgalgrid - deflelem[:, 1]
                        lgalprim = gdat.lgalgrid - deflelem[:, 0]
                        dictelem[l]['relm'][k] = np.mean(abs(cntp['lens'][0, :, 0] - cntplensobjt(bgalprim, lgalprim, grid=False).flatten()))
                        
                        
                        dictelem[l]['relk'][k] = dictelem[l]['relm'][k] / dictelem[l]['defs'][k] * gdat.sizepixl
                        dictelem[l]['reln'][k] = dictelem[l]['rele'][k] / dictelem[l]['defs'][k] * gdat.sizepixl
                        dictelem[l]['reld'][k] = retr_rele(gdat, gdat.cntpdata[0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], \
                                                                                                         dictelem[l]['defs'][k], asca, acut, gdat.indxpixl)
                        dictelem[l]['relc'][k] = retr_rele(gdat, cntp['lens'][0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], \
                                                       dictelem[l]['defs'][k], asca, acut, gdat.indxpixl, absv=False) / dictelem[l]['defs'][k] * gdat.sizepixl
               
        ### distribution of element parameters and features
        #### calculate the model filter
        listindxelemfilt = [[[] for l in indxpopl] for namefilt in gdat.listnamefilt]
        for k, namefilt in enumerate(gdat.listnamefilt):
            for l in indxpopl:
                if namefilt == '':
                    listindxelemfilt[k][l] = np.arange(numbelem[l])
                if namefilt == 'imagbndr':
                    listindxelemfilt[k][l] = np.where((np.fabs(dictelem[l]['lgal']) < gdat.maxmgangdata) & (np.fabs(dictelem[l]['bgal']) < gdat.maxmgangdata))[0]
                if namefilt == 'deltllik':
                    listindxelemfilt[k][l] = np.where(dictelem[l]['deltllik'] > 0.5 * numbcomp[l])[0]
                if namefilt == 'nrel':
                    listindxelemfilt[k][l] = np.where(dictelem[l]['reln'] > 0.3)[0]
    
        listnamefeatsele = getattr(gdat, strgmodl + 'listnamefeatsele')
        
        # histograms of element features
        for l in indxpopl:
            #### one dimensional
            for strgfeat in liststrgfeat[l]:
                if strgfeat == 'spec':
                    temp = np.zeros((gdat.numbbinsplot, gdat.numbener))
                else:
                    temp = np.zeros(gdat.numbbinsplot)
               
                if strgfeat[:-4] == 'etag':
                    continue
                dictelem[l]['hist' + strgfeat] = temp
                if strgfeat == 'specplot' or strgfeat == 'deflprof':
                    continue
                elif strgfeat == 'spec':
                    for i in gdat.indxener:
                        dictelem[l]['hist' + strgfeat][:, i] = np.histogram(dictelem[l]['spec'][i, listindxelemfilt[0][l]], gdat.binsspec)[0]
                elif strgfeat == 'cnts':
                    dictelem[l]['hist' + strgfeat] = np.histogram(dictelem[l]['cnts'][listindxelemfilt[0][l]], gdat.binscnts)[0]
                elif not (strgfeat == 'curv' and spectype[l] != 'curv' or strgfeat == 'expc' \
                                                        and spectype[l] != 'expc' or strgfeat.startswith('sindarry') and \
                                                                                                                         spectype[l] != 'colr'):
                    bins = getattr(gdat, 'bins' + strgfeat)
                    if len(dictelem[l][strgfeat]) > 0 and len(listindxelemfilt[0][l]) > 0:
                        dictelem[l]['hist' + strgfeat] = np.histogram(dictelem[l][strgfeat][listindxelemfilt[0][l]], bins)[0]

            #### two dimensional
            dictelemtdim = dict()
            for l0 in indxpopl:
                for strgfrst in liststrgfeat[l0]:
                    for strgseco in liststrgfeat[l0]:
                        
                        if strgfrst == 'spec' or strgfrst == 'specplot' or strgfrst == 'deflprof' or \
                                strgseco == 'spec' or strgseco == 'specplot' or strgseco == 'deflprof':
                            continue
                        
                        if not checstrgfeat(strgfrst, strgseco):
                            continue

                        strgtemp = 'hist' + strgfrst + strgseco + 'pop%d' % l0
                        dictelemtdim[strgtemp] = np.zeros((gdat.numbbinsplot, gdat.numbbinsplot))
            
                        binsfrst = getattr(gdat, 'bins' + strgfrst)
                        binsseco = getattr(gdat, 'bins' + strgseco)
                        if len(dictelem[l0][strgfrst]) > 0 and len(dictelem[l0][strgseco]) > 0:
                            dictelemtdim[strgtemp] = np.histogram2d(dictelem[l0][strgfrst][listindxelemfilt[0][l0]], \
                                                                    dictelem[l0][strgseco][listindxelemfilt[0][l0]], [binsfrst, binsseco])[0]
                        strg = strgpfix + strgtemp
                        
                        setattr(gdatobjt, strg, dictelemtdim[strgtemp])
            
            ### priors on element parameters and features
            dictelem[l]['hist' + strgfeat + 'prio'] = np.empty(gdat.numbbinsplotprio)
            for strgfeat in liststrgfeatprio[l]:
                strgpdfn = liststrgpdfnprio[l][liststrgfeatprio[l].index(strgfeat)]
                minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                bins = getattr(gdat, 'bins' + strgfeat)
                delt = getattr(gdat, 'delt' + strgfeat)
                
                xdatplot = getattr(gdat, 'mean' + strgfeat)
                xdat = getattr(gdat, strgmodl + 'mean' + strgfeat + 'prio')
                deltprio = getattr(gdat, strgmodl + 'delt' + strgfeat + 'prio')
                
                booltemp = False
                if strgpdfn.startswith('expo') or strgpdfn.startswith('dexp'):
                    if strgpdfn.startswith('expo'):
                        if strgpdfn == 'expo':
                            sexp = getattr(gdat, strgmodl + 'gangdistsexppop%d' % l)
                        else:
                            sexp = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distscal')[l]]
                        pdfn = pdfn_expo(xdat, maxm, sexp)
                    if strgpdfn.startswith('dexp'):
                        pdfn = pdfn_dnp.exp(xdat, maxm, scal)
                    booltemp = True
                if strgpdfn.startswith('self') or strgpdfn.startswith('logt'):
                    minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                    maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                    if strgpdfn.startswith('self'):
                        pdfn = 1. / (maxm - minm) + np.zeros_like(xdat)
                    else:
                        pdfn = 1. / (np.log(maxm) - np.log(minm)) + np.zeros_like(xdat)
                    booltemp = True
                # temp 
                if strgpdfn.startswith('powrslop'):
                    slop = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
                    pdfn = pdfn_powr(xdat, minm, maxm, slop)
                    booltemp = True
                if strgpdfn.startswith('dpowslopbrek'):
                    brek = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distbrek')[l]]
                    sloplowr = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distsloplowr')[l]]
                    slopuppr = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslopuppr')[l]]
                    pdfn = pdfn_dpow(xdat, minm, maxm, brek, sloplowr, slopuppr)
                    booltemp = True
                if strgpdfn == 'lnormeanstdv':
                    meanlnor = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distmean')[l]]
                    stdvlnor = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'diststdv')[l]]
                    pdfn = pdfn_lnor(xdat, meanlnor, stdvlnor)
                    booltemp = True

                if strgpdfn.startswith('igam'):
                    cutf = getattr(gdat, 'cutf' + strgfeat)
                    pdfn = pdfn_igam(xdat, slop, cutf)
                    booltemp = True
                if strgpdfn.startswith('gaus'):
                    # this does not work for mismodeling
                    meanvarb = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distmean')[l]]
                    stdv = samp[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'diststdv')[l]]
                    if strgfeat == 'expc' and spectype[l] == 'expc':
                        pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                    else:
                        pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                    booltemp = True
                
                if booltemp:
                    dictelem[l]['hist' + strgfeat + 'prio'] = meanelem[l] * pdfn * np.interp(xdat, xdatplot, delt)
                
                setattr(gdatobjt, strgpfix + 'hist' + strgfeat + 'pop%dprio' % l, dictelem[l]['hist' + strgfeat + 'prio'])
                if strgmodl == 'true':
                    setattr(gdatobjt, 'refrhist' + strgfeat + 'pop%dprio' % l, dictelem[l]['hist' + strgfeat + 'prio'])
    
    if numbtrap > 0:
        for l in indxpopl:
            if elemtype[l] == 'lens':
                if numbelempopl[l] > 0:
                    ## total truncated mass of the subhalo as a cross check
                    # temp -- generalize
                    asca = dictelem[l]['asca']
                    acut = dictelem[l]['acut']
                    factmcutfromdefs = retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut) 
                    masssubh = np.array([np.sum(factmcutfromdefs * dictelem[l]['defs'])])
    
    ## derived variables as a function of other derived variables
    if numbtrap > 0:
        for l in indxpopl:
            if elemtype[l].startswith('lghtpntspuls'):
                massshel = np.empty(gdat.numbanglhalf)
                for k in gdat.indxanglhalf:
                    indxelemshel = np.where((gdat.binsanglhalf[k] < dictelem[l]['gang']) & (dictelem[l]['gang'] < gdat.binsanglhalf[k+1]))
                    massshel[k] = np.sum(dictelem[l]['mass'][indxelemshel])
                setattr(gdatobjt, strgpfix + 'massshelpop%d' % l, massshel)
            
    if lensmodltype != 'none' or numbtrap > 0 and 'lens' in elemtype:
        # find the host, subhalo masses and subhalo mass fraction as a function of halo-centric radius
        mdencrit = getattr(gdat, strgmodl + 'mdencrit')
        adishost = getattr(gdat, strgmodl + 'adishost')
        adishostsour = getattr(gdat, strgmodl + 'adishostsour')
        massfrombein = getattr(gdat, strgmodl + 'massfrombein')
        listnametemp = ['delt', 'intg']
        listnamevarbmass = []
        listnamevarbmassscal = []
        listnamevarbmassvect = []
        for e in indxsersfgrd:
            if lensmodltype == 'host' or lensmodltype == 'full':
                listnamevarbmassscal += ['masshosttotl']
                for strgtemp in listnametemp:
                    listnamevarbmassvect.append('masshostisf%d' % e + strgtemp)
                    listnamevarbmassscal.append('masshostisf%d' % e + strgtemp + 'bein')
        if numbtrap > 0 and 'lens' in elemtype:
            listnamevarbmassscal.append('masssubhtotl')
            listnamevarbmassscal.append('fracsubhtotl')
            for strgtemp in listnametemp:
                listnamevarbmassvect.append('masssubh' + strgtemp)
                listnamevarbmassvect.append('fracsubh' + strgtemp)
                listnamevarbmassscal.append('masssubh' + strgtemp + 'bein')
                listnamevarbmassscal.append('fracsubh' + strgtemp + 'bein')

        for name in listnamevarbmassvect:
            dicttert[name] = np.zeros(gdat.numbanglhalf)
            if 'isf' in name:
                indxisfrtemp = int(name.split('isf')[1][0])
            angl = np.sqrt((gdat.meanlgalcartmesh - lgalhost[indxisfrtemp])**2 + (gdat.meanbgalcartmesh - bgalhost[indxisfrtemp])**2).flatten()
            for k in gdat.indxanglhalf:
                if name[4:8] == 'host':
                    convtemp = conv[:]
                if name[4:8] == 'subh':
                    convtemp = convelem[:]
                
                if name.endswith('delt'):
                    indxpixl = np.where((gdat.binsanglhalf[k] < angl) & (angl < gdat.binsanglhalf[k+1]))[0]
                    dicttert[name][k] = 1e6 * np.sum(convtemp[indxpixl]) * mdencrit * \
                                                gdat.apix * adishost**2 / 2. / np.pi * gdat.deltanglhalf[k] / gdat.meananglhalf[k]
                if name.endswith('intg'):
                    indxpixl = np.where(angl < gdat.meananglhalf[k])[0]
                    dicttert[name][k] = np.sum(convtemp[indxpixl]) * mdencrit * gdat.apix * adishost**2
                
                if name[:4] == 'frac':
                    masshosttotl = 0.
                    for e in indxsersfgrd:
                        masshosttotl += dicttert['masshostisf%d' % e + name[-4:]][k]
                    if masshosttotl != 0.:
                        dicttert['fracsubh' + name[8:]][k] = dicttert['masssubh' + name[8:]][k] / masshosttotl
            setattr(gdatobjt, strgpfix + name, dicttert[name])
            
            # interpolate the host, subhalo masses and subhalo mass fraction at the Einstein radius and save it as a scalar variable
            dicttert[name + 'bein'] = np.interp(beinhost, gdat.meananglhalf, dicttert[name])
            setattr(gdatobjt, strgpfix + name + 'bein', dicttert[name + 'bein'])
        
    if numbtrap > 0:
        ## copy element features to the global object
        feat = [[] for l in indxpopl]
        for l in indxpopl:
            feat[l] = dict()
            for strgfeat in liststrgfeat[l]:
                if strgfeat[:-4] == 'etag':
                    continue
                if len(dictelem[l][strgfeat]) > 0:
                    if strgmodl == 'true':
                        shap = list(np.ones(dictelem[l][strgfeat].ndim, dtype=int))
                        feat[l][strgfeat] = np.tile(dictelem[l][strgfeat], [3] + shap)
                    if strgmodl == 'fitt':
                        feat[l][strgfeat] = dictelem[l][strgfeat]
                setattr(gdatobjt, strgpfix + 'hist' + strgfeat + 'pop%d' % l, dictelem[l]['hist' + strgfeat])
                    
        for strgfeat in liststrgfeattotl:
            feattemp = [[] for l in indxpopl]
            for l in indxpopl:
                if strgfeat in liststrgfeat[l]:
                    if strgfeat in feat[l]:
                        feattemp[l] = feat[l][strgfeat]
                    else:
                        feattemp[l] = np.array([])
            setattr(gdatobjt, strgpfix + strgfeat, feattemp)
        
    # copy true state to the reference state
    if strgmodl == 'true':
        for name, valu in deepcopy(gdat.__dict__).iteritems():
            if name.startswith('true'):
                #indx = name.find('pop')
                #if indx != -1 and not name.endswith('pop') and name[indx+3].isdigit():
                #    namerefr = name.replace('pop%s' % name[indx+3], 'ref%s' % name[indx+3])
                #else:
                #    namerefr = name
                #namerefr = name
                #namerefr = namerefr.replace('true', 'refr')
                name = name.replace('true', 'refr')
                setattr(gdat, name, valu)
    
    if numbtrap > 0 and gdat.priofactdoff != 0.:
        if strgmodl == 'true':
            for q in gdat.indxrefr:
                for strgfeat in gdat.refrliststrgfeat[q]:
                    if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':
                        continue
                    refrhist = getattr(gdat, 'truehist' + strgfeat + 'pop%d' % q)
                    
                    indxelempars = np.where(dictelem[q]['deltllik'] > 2.5)[0]
                    
                    ptfn = np.zeros(gdat.numbbinsplot) - 1.
                    refrhistpars = np.zeros(gdat.numbbinsplot) - 1.
                    
                    indxrefrgood = np.where(refrhist > 0)[0]
                    ptfn[indxrefrgood] = 0.
                    refrhistpars[indxrefrgood] = 0.
                    
                    bins = getattr(gdat, 'bins' + strgfeat)
                    if len(indxelempars) > 0:
                        refrhistpars = np.histogram(dictelem[q][strgfeat][indxelempars], bins=bins)[0].astype(float)
                        if indxrefrgood.size > 0:
                            ptfn[indxrefrgood] = refrhistpars[indxrefrgood] / refrhist[indxrefrgood]
                    
                    setattr(gdatobjt, strgpfix + 'histpars' + strgfeat + 'pop%d' % q, refrhistpars)
                    setattr(gdatobjt, strgpfix + 'ptfn' + strgfeat + 'pop%d' % q, ptfn)
        
        if gdat.rtagmock is not None and gdat.datatype == 'inpt' or gdat.datatype == 'mock':
            if gdat.truenumbtrap > 0:
                for l in indxpopl:
                    for strgfeat in liststrgfeat[l]:
                        if strgfeat == 'spec' or strgfeat == 'specplot' or strgfeat == 'deflprof':# or strgfeat.startswith('aerr'):
                            continue
                        if strgfeat in gdat.trueliststrgfeat[l]:
                            hist = getattr(gdatobjt, strgpfix + 'hist' + strgfeat + 'pop%d' % l)
                            ptfn = getattr(gdat, 'trueptfn' + strgfeat + 'pop%d' % l)
                            histptfn = hist / ptfn
                            setattr(gdatobjt, strgpfix + 'histptfn' + strgfeat + 'pop%d' % l, histptfn)

    ### Exculusive comparison with the true state
    if strgmodl == 'fitt' and gdat.datatype == 'mock':
        if lensmodltype != 'none':
            numbsingcomm = min(deflsing.shape[2], gdat.truedeflsing.shape[2])
            deflsingresi = deflsing[0, ..., :numbsingcomm] - gdat.truedeflsing[..., :numbsingcomm]
            deflsingresimgtd = np.sqrt(np.sum(deflsingresi**2, axis=1))
            deflsingresiperc = 100. * deflsingresimgtd / gdat.truedeflsingmgtd[..., :numbsingcomm]
            setattr(gdatobjt, strgpfix + 'numbsingcomm', numbsingcomm)
            setattr(gdatobjt, strgpfix + 'deflsingresi', deflsingresi)
            truedeflmgtd = getattr(gdat, 'truedeflmgtd')
            truedefl = getattr(gdat, 'truedefl')
            deflresi = defl - truedefl
            deflresimgtd = np.sqrt(np.sum(deflresi**2, axis=1))
            deflresiperc = 100. * deflresimgtd / truedeflmgtd
            setattr(gdatobjt, strgpfix + 'deflresi', deflresi)
            setattr(gdatobjt, strgpfix + 'deflresimgtd', deflresimgtd)
            if numbtrap > 0:
                trueconvelem = getattr(gdat, 'trueconvelem')
                convelemresi = convelem[:] - trueconvelem
                convelemresiperc = 100. * convelemresi / trueconvelem
                setattr(gdatobjt, strgpfix + 'convelemresi', convelemresi)
                setattr(gdatobjt, strgpfix + 'convelemresiperc', convelemresiperc)
            truemagn = getattr(gdat, 'truemagn')
            magnresi = magn[:] - truemagn
            magnresiperc = 100. * magnresi / truemagn
            setattr(gdatobjt, strgpfix + 'magnresi', magnresi)
            setattr(gdatobjt, strgpfix + 'magnresiperc', magnresiperc)
    
    if numbtrap > 0:
        # correlate the catalog sample with the reference catalog
        if gdat.refrinfo and not (strgmodl == 'true' and gdat.datatype == 'mock') and gdat.boolasscrefr:
            for q in gdat.indxrefr:
                for l in gdat.fittindxpopl:
                    if gdat.refrnumbelem[q] > 0:
                        cmpl = np.array([float(len(indxelemrefrasschits[q][l])) / gdat.refrnumbelem[q]])
                        if gdat.booldiagmode:
                            if cmpl > 1. or cmpl < 0.:
                                raise Exception('')
                    else:
                        cmpl = np.array([-1.])
                    setattr(gdatobjt, strgpfix + 'cmplpop%dpop%d' % (l, q), cmpl)
                    if numbelem[l] > 0:
                        fdis = np.array([float(indxelemfittasscfals[q][l].size) / numbelem[l]])
                        if gdat.booldiagmode:
                            if fdis > 1. or fdis < 0.:
                                raise Exception('')
                    else:
                        fdis = np.array([-1.])
                    setattr(gdatobjt, strgpfix + 'fdispop%dpop%d' % (q, l), fdis)
                        
            # collect the associated fitting element feature for each reference element
            featrefrassc = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
            for q in gdat.indxrefr:
                for l in gdat.fittindxpopl:
                    featrefrassc[q][l] = dict()
                    for strgfeat in gdat.refrliststrgfeat[q]:
                        if not strgfeat in liststrgfeat[l] or strgfeat in gdat.refrliststrgfeatonly[q][l]:
                            continue
                        if isinstance(dictelem[l][strgfeat], np.ndarray) and dictelem[l][strgfeat].ndim > 1:
                            continue
                        featrefrassc[q][l][strgfeat] = np.zeros(gdat.refrnumbelem[q]) + np.nan
                        if len(indxelemrefrasschits[q][l]) > 0 and len(dictelem[l][strgfeat]) > 0:
                            featrefrassc[q][l][strgfeat][indxelemrefrasschits[q][l]] = dictelem[l][strgfeat][indxelemfittasschits[q][l]]
                        name = strgpfix + strgfeat + 'asscpop%dpop%d' % (q, l)
                        setattr(gdatobjt, name, featrefrassc[q][l][strgfeat])
            
            # completeness
            for l in gdat.fittindxpopl:
                for q0 in gdat.indxrefr:
                    if gdat.refrnumbelem[q0] == 0:
                        continue
                
                    for (strgfeatfrst, strgfeatfrsttagg) in zip(gdat.refrliststrgfeat[q0], gdat.refrliststrgfeattagg[q0][l]):
                        
                        if strgfeatfrst.startswith('etag'):
                            continue
                                
                        if strgfeatfrst == 'spec' or strgfeatfrst == 'specplot':
                            continue
                        
                        refrfeatfrst = getattr(gdat, 'refr' + strgfeatfrst)
                        binsfeatfrst = getattr(gdat, 'bins' + strgfeatfrst)
                            
                        for (strgfeatseco, strgfeatsecotagg) in zip(gdat.refrliststrgfeat[q0], gdat.refrliststrgfeattagg[q0][l]):
                            if strgfeatfrst == strgfeatseco:
                                continue
                            
                            if strgfeatseco.startswith('etag'):
                                continue
                            
                            if strgfeatseco == 'spec' or strgfeatseco == 'specplot':
                                continue
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                            
                            # temp -- the size of the cmpl np.array should depend on strgmodl
                            cmpltdim = np.zeros((gdat.numbbinsplot, gdat.numbbinsplot)) - 1.
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                            
                            if len(indxelemrefrasschits[q0][l]) > 0:
                                strgfeattdim = strgfeatfrst + strgfeatseco + 'pop%d' % q0
                                refrhistfeattdim = getattr(gdat, 'refrhist' + strgfeattdim)
                                refrfeatseco = getattr(gdat, 'refr' + strgfeatseco)
                                binsfeatseco = getattr(gdat, 'bins' + strgfeatseco)
                                
                                refrhistfeattdimassc = np.histogram2d(refrfeatfrst[q0][0, indxelemrefrasschits[q0][l]], \
                                                                   refrfeatseco[q0][0, indxelemrefrasschits[q0][l]], bins=(binsfeatfrst, binsfeatseco))[0]
                                indxgood = np.where(refrhistfeattdim != 0.)
                                if indxgood[0].size > 0:
                                    cmpltdim[indxgood] = refrhistfeattdimassc[indxgood].astype(float) / refrhistfeattdim[indxgood]
                                    if gdat.booldiagmode:
                                        if np.where((cmpltdim[indxgood] > 1.) | (cmpltdim[indxgood] < 0.))[0].size > 0:
                                            print 'Warning! Completeness went outside [0, 1]'
                                            print 'strgfeatfrst'
                                            print strgfeatfrst
                                            print 'strgfeatseco'
                                            print strgfeatseco
                                            print 'binsfeatfrst'
                                            print binsfeatfrst
                                            print 'binsfeatseco'
                                            print binsfeatseco
                                            print 'refrfeatfrst[q0]'
                                            print refrfeatfrst[q0]
                                            print 'refrfeatseco[q0]'
                                            print refrfeatseco[q0]
                                            print 'indxelemrefrasschits[q0]'
                                            print indxelemrefrasschits[q0]
                                            print 'strgfeattdim'
                                            print strgfeattdim
                                            print 'refrhistfeattdim'
                                            print refrhistfeattdim
                                            summgene(refrhistfeattdim)
                                            print 'refrhistfeattdimassc'
                                            print refrhistfeattdimassc
                                            summgene(refrhistfeattdimassc)
                                            print 'indxgood'
                                            print indxgood
                                            print
                                            raise Exception('')
                            
                            strg = strgpfix + 'cmpl' + strgfeatfrst + strgfeatseco + 'pop%dpop%d' % (l, q0)
                            setattr(gdatobjt, strg, cmpltdim)

                        cmplfeatfrst = np.zeros(gdat.numbbinsplot) - 1.
                        if len(indxelemrefrasschits[q0][l]) > 0:
                            refrhistfeatfrst = getattr(gdat, 'refrhist' + strgfeatfrst + 'pop%d' % q0)
                            binsfeatfrst = getattr(gdat, 'bins' + strgfeatfrst)
                            refrhistfeatfrstassc = np.histogram(refrfeatfrst[q0][0, indxelemrefrasschits[q0][l]], bins=binsfeatfrst)[0]
                            indxgood = np.where(refrhistfeatfrst != 0.)[0]
                            if indxgood.size > 0:
                                cmplfeatfrst[indxgood] = refrhistfeatfrstassc[indxgood].astype(float) / refrhistfeatfrst[indxgood]
                                if gdat.booldiagmode:
                                    if np.where((cmplfeatfrst[indxgood] > 1.) | (cmplfeatfrst[indxgood] < 0.))[0].size > 0:
                                        raise Exception('')
                       
                        setattr(gdatobjt, strgpfix + 'cmpl' + strgfeatfrst + 'pop%dpop%d' % (l, q0), cmplfeatfrst)
                        
            # false discovery rate
            for q in gdat.indxrefr:
                for l0 in gdat.fittindxpopl:
                    for strgfeatfrst in liststrgfeatodim[l0]:
                        
                        #if (not strgfeatfrst in gdat.refrliststrgfeat[q]) and (not strgfeatfrst in gdat.refrliststrgfeattagg[q][l0]):
                        #if not strgfeatfrst in gdat.refrliststrgfeattagg[q][l0]:
                        #    continue
                
                        if strgfeatfrst.startswith('etag') or strgfeatfrst.startswith('aerr'):
                            continue
                                
                        if strgfeatfrst == 'spec' or strgfeatfrst == 'specplot' or strgfeatfrst == 'deflprof':
                            continue
                        
                        binsfeatfrst = getattr(gdat, 'bins' + strgfeatfrst)
                        for strgfeatseco in liststrgfeatodim[l0]:
                            
                            #if (not strgfeatseco in gdat.refrliststrgfeat[q]) and (not strgfeatseco in gdat.refrliststrgfeattagg[q][l0]):
                            #if not strgfeatseco in gdat.refrliststrgfeattagg[q][l0]:
                            #    continue
    
                            if strgfeatseco.startswith('aerr'):
                                continue
                            
                            strgfeattdim = strgfeatfrst + strgfeatseco + 'pop%d' % l0
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                            
                            # temp -- the size of the fdis np.array should depend on strgmodl
                            fdistdim = np.zeros((gdat.numbbinsplot, gdat.numbbinsplot))
                            
                            if len(indxelemrefrasschits[q][l0]) > 0 and len(dictelem[l0][strgfeatseco]) > 0 and len(dictelem[l0][strgfeatfrst]) > 0: 
                                fitthistfeattdim =  dictelemtdim['hist' + strgfeattdim]
                                binsfeatseco = getattr(gdat, 'bins' + strgfeatseco)
                                
                                fitthistfeattdimfals = np.histogram2d(dictelem[l0][strgfeatfrst][indxelemfittasscfals[q][l0]], \
                                                      dictelem[l0][strgfeatseco][indxelemfittasscfals[q][l0]], bins=(binsfeatfrst, binsfeatseco))[0]
                                indxgood = np.where(fitthistfeattdim != 0.)
                                if indxgood[0].size > 0:
                                    fdistdim[indxgood] = fitthistfeattdimfals[indxgood].astype(float) / fitthistfeattdim[indxgood]
                                    if gdat.booldiagmode:
                                        if np.where((fdistdim[indxgood] > 1.) | (fdistdim[indxgood] < 0.))[0].size > 0:
                                            raise Exception('')
                            
                            setattr(gdatobjt, strgpfix + 'fdis' + strgfeatfrst + strgfeatseco + 'pop%dpop%d' % (q, l0), fdistdim)
                    
                        fdisfeatfrst = np.zeros(gdat.numbbinsplot)
                        if len(indxelemrefrasschits[q][l0]) > 0 and len(dictelem[l0][strgfeatfrst]) > 0:
                            binsfeatfrst = getattr(gdat, 'bins' + strgfeatfrst)
                            fitthistfeatfrstfals = np.histogram(dictelem[l0][strgfeatfrst][indxelemfittasscfals[q][l0]], bins=binsfeatfrst)[0]
                            fitthistfeatfrst = getattr(gdatobjt, strgpfix + 'hist' + strgfeatfrst + 'pop%d' % l0)
                            indxgood = np.where(fitthistfeatfrst != 0.)[0]
                            if indxgood.size > 0:
                                fdisfeatfrst[indxgood] = fitthistfeatfrstfals[indxgood].astype(float) / fitthistfeatfrst[indxgood]
                                if gdat.booldiagmode:
                                    if np.where((fdisfeatfrst[indxgood] > 1.) | (fdisfeatfrst[indxgood] < 0.))[0].size > 0:
                                        
                                        print 'Warning! FDR went outside [0, 1]'
                                        print 'strgfeatfrst'
                                        print strgfeatfrst
                                        #print 'strgfeatseco'
                                        #print strgfeatseco
                                        print 'binsfeatfrst'
                                        print binsfeatfrst
                                        print 'fitthistfeatfrst'
                                        print fitthistfeatfrst
                                        print 'fitthistfeatfrstfals'
                                        print fitthistfeatfrstfals
                                        print
                                        raise Exception('')
                        
                        setattr(gdatobjt, strgpfix + 'fdis' + strgfeatfrst + 'pop%dpop%d' % (q, l0), fdisfeatfrst)
    
    
        namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
      
        # temp
        if strgmodl == 'true' and gdat.verbtype > 0:
            for l in indxpopl:
                for strgfeat in liststrgfeat[l]:
                    #minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                    #maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                    minm = getattr(gdat, 'minm' + strgfeat)
                    maxm = getattr(gdat, 'maxm' + strgfeat)
                    if np.where(minm > dictelem[l][strgfeat])[0].size > 0 or np.where(maxm < dictelem[l][strgfeat])[0].size > 0:
                        print 'Warning: element feature outside the plot limits.'
                        print 'l'
                        print l
                        print 'Feature: '
                        print strgfeat
                        print 'Plot minmimum'
                        print minm
                        print 'Plot maxmimum'
                        print maxm
                        if len(dictelem[l][strgfeat]) > 0:
                            print 'Feature np.minimum'
                            print np.amin(dictelem[l][strgfeat])
                            print 'Feature np.maximum'
                            print np.amax(dictelem[l][strgfeat])
                            print
                        if strgfeat == namefeatampl[l] and strgfeat in liststrgcomp[l]:
                            indxcomptemp = liststrgcomp[l].index(strgfeat)
                            if (listscalcomp[l][indxcomptemp] != 'gaus' and not listscalcomp[l][indxcomptemp].startswith('lnor')):
                                raise Exception('')
    stopchro(gdat, gdatmodi, 'tert')
    
        
def retr_lprielem(gdat, strgmodl, l, g, strgfeat, strgpdfn, samp, dictelem, numbelem):
    
    if strgpdfn == 'self':
        minmfeat = getattr(gdat, strgmodl + 'minm' + strgfeat)
        maxmfeat = getattr(gdat, strgmodl + 'maxm' + strgfeat)
        lpri = numbelem[l] * np.log(1. / (maxmfeat - minmfeat))
    if strgpdfn == 'logt':
        lpri = retr_lprilogtdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, samp, l)
    if strgpdfn == 'gaus':
        lpri = retr_lprigausdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, samp, l)
    if strgpdfn == 'dexpscal':
        maxmbgal = getattr(gdat, strgmodl + 'maxmbgal')
        indxfixpbgaldistscal = getattr(gdat, strgmodl + 'indxfixpbgaldistscalpop%d' % l)
        lpri = np.sum(np.log(pdfn_dnp.exp(dictelem[l]['bgal'], maxmbgal, samp[indxfixpbgaldistscal]))) 
    if strgpdfn == 'exposcal':
        maxmgang = getattr(gdat, strgmodl + 'maxmgang')
        gang = retr_gang(dictelem[l]['lgal'], dictelem[l]['bgal'])
        indxfixpgangdistscal = getattr(gdat, strgmodl + 'indxfixpgangdistscalpop%d' % l)
        lpri = np.sum(np.log(pdfn_expo(gang, maxmgang, samp[indxfixpgangdistscal]))) 
        lpri = -numbelem[l] * np.log(2. * pi) 
    if strgpdfn == 'tmpl':
        lpri = np.sum(lpdfspatprioobjt(dictelem[l]['bgal'], dictelem[l]['lgal'], grid=False))
    if strgpdfn == 'powrslop':
        lpri = retr_lpripowrdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, samp, l)
    if strgpdfn == 'dpowslopbrek':
        lpri = retr_lpridpowdist(gdat, strgmodl, dictelem[l][strgfeat], strgfeat, samp, l)
    if strgpdfn == 'dsrcexpo':
        lpri += -np.sum(np.sqrt((dictelem[l]['lgal'] - lgalsour)**2 + (dictelem[l]['bgal'] - bgalsour)**2) / \
                                                                getattr(gdat, strgmodl + 'dsrcdistsexppop%d' % l))
    if strgpdfn == 'tmpl':
        if strgpdfn.endswith('cons'):
            pdfnspatpriotemp = getattr(gdat, strgmodl + 'pdfnspatpriotemp')
            spatdistcons = samp[getattr(gdat, strgmodl + 'indxfixpspatdistcons')]
            lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons)
            lpdfspatpriointp = lpdfspatprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
            lpdfspatpriointp = lpdfspatpriointp.T
            setattr(gdatobjt, strgpfix + 'lpdfspatpriointp', lpdfspatpriointp)
            setattr(gdatobjt, strgpfix + 'lpdfspatprioobjt', lpdfspatprioobjt)
        else:
            lpdfspatprioobjt = gdat.fittlpdfspatprioobjt
    
    return lpri


def checstrgfeat(strgfrst, strgseco):

    numbfrst = len(strgfrst)
    numbseco = len(strgseco)
    numb = min(numbfrst, numbseco)
    if strgfrst[:numb] < strgseco[:numb]:
        booltemp = True
    elif strgfrst[:numb] == strgseco[:numb]:
        if numbfrst >= numbseco:
            booltemp = False
        else:
            booltemp = True
    else:
        booltemp = False

    return booltemp


def retr_pathoutprtag(rtag):
    
    pathoutprtag = os.environ["PCAT_DATA_PATH"] + '/data/outp/' + rtag + '/'
    
    return pathoutprtag


def proc_finl(gdat=None, rtag=None, strgpdfn='post', listnamevarbproc=None, forcplot=False):
    
    gdatmock = None
    
    print 'proc_finl()'

    if rtag is None:
        rtag = gdat.rtag
    
    # determine if the final-processing if nominal or tiling
    if isinstance(rtag, list):
        listrtagmodi = rtag
        rtagfinl = tdpy.util.retr_strgtimestmp() + rtag[0][15:] + 'tile'
        booltile = True
    else:
        listrtagmodi = [rtag]
        rtagfinl = rtag
        booltile = False
    
    # determine of the gdatfinl object is available 
    boolgdatfinl = chec_statfile(rtagfinl, 'gdatfinl', strgpdfn)
    boolgdatfinlgood = False
    if boolgdatfinl:
        print 'Final-processing has been performed previously.'
        pathoutprtag = retr_pathoutprtag(rtagfinl)
        path = pathoutprtag + 'gdatfinl' + strgpdfn
        try:
            gdat = readfile(path) 
            boolgdatfinlgood = True
        except:
            print 'gdatfinl object is corrupted.'

    if boolgdatfinl and boolgdatfinlgood:
        # read gdatfinl
        pathoutprtag = retr_pathoutprtag(rtagfinl)
        path = pathoutprtag + 'gdatfinl' + strgpdfn
        gdatfinl = readfile(path) 
        
        if gdatfinl.fittnumbtrap > 0:
            if gdatfinl.datatype == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:
                    if gdatfinl.rtagmock is not None:
                        path = gdatfinl.pathoutprtagmock + 'gdatfinlpost'
                        gdatmock = readfile(path)
                    
    else:
        
        if booltile:
            gdatfinltile = tdpy.util.gdatstrt()
        
        indxrtaggood = []
        liststrgtile = []
        listrtaggood = []
        indxtiletemp = 0
        for n, rtagmodi in enumerate(listrtagmodi):
            
            # read gdatinit
            boolgdatinit = chec_statfile(rtagmodi, 'gdatinit', '')
            if not boolgdatinit:
                if booltile:
                    print 'Initial global object not found. Skipping...'
                    continue
                else:
                    print 'Initial global object not found. Quitting...'
                    return
            
            pathoutprtag = retr_pathoutprtag(rtagmodi)
            path = pathoutprtag + 'gdatinit'
            
            gdatinit = readfile(path) 
            if booltile:
                gdatfinltile = gdatinit
                gdatfinl = gdatinit
            else:
                gdatfinl = gdatinit

            if gdatinit.mockonly:
                print 'Mock only run. Quitting final-processing...'
                return

            # read gdatmodi
            boolgdatmodi = chec_statfile(rtagmodi, 'gdatmodi', strgpdfn)
            if not boolgdatmodi:
                print 'Modified global object not found. Quitting final-processing...'
                return
        
            ## list of other parameters to be flattened
            gdatinit.liststrgvarbarryflat = deepcopy(gdatinit.liststrgvarbarry)
            for strg in ['memoresi']:
                gdatinit.liststrgvarbarryflat.remove(strg)
   
            listsamp = np.empty((gdatinit.numbsamptotl, gdatinit.fittnumbpara))
            
            if booltile:
                gdatfinltile.pathoutprtag = retr_pathoutprtag(rtagfinl)
                numbsamptotlrsmp = gdatinit.numbsamptotl
                indxsamptotlrsmp = choice(gdatinit.indxsamptotl, size=gdatinit.numbsamptotl, replace=False)
            
            pathoutprtagmodi = retr_pathoutprtag(rtagmodi)
            
            if rtagmodi.split('_')[-2][-4:] in liststrgtile:
                print 'Tile previously processed. Skipping...'
                print
                continue

            # aggregate samples from the chains
            if gdatinit.verbtype > 0:
                print 'Reading gdatmodi objects from all processes...'
                timeinit = gdatinit.functime()
            
            listgdatmodi = []
            for k in gdatinit.indxproc:
                path = pathoutprtagmodi + 'gdatmodi%04d' % k + strgpdfn
                listgdatmodi.append(readfile(path))
            
            if gdatinit.verbtype > 0:
                timefinl = gdatinit.functime()
                print 'Done in %.3g seconds.' % (timefinl - timeinit)
            
            if gdatinit.fittnumbtrap > 0:
                if len(getattr(listgdatmodi[0], 'list' + strgpdfn + 'indxelemfull')) == 0:
                    print 'Found an empty element list. Skipping...'
                    print
                    continue
            
            if gdatinit.verbtype > 0:
                print 'Accumulating np.arrays...'
                timeinit = gdatinit.functime()
            
            for strgvarb in gdatinit.liststrgvarbarry:
                for k in gdatinit.indxproc:
                    if k == 0:
                        shap = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb).shape
                        shap = [shap[0], gdatinit.numbproc] + list(shap[1:])
                        temp = np.zeros(shap) - 1
                    if len(shap) > 2:
                        temp[:, k, :] = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb)
                    else:
                        temp[:, k] = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb)
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, temp)
            
            if gdatfinl.verbtype > 0:
                timefinl = gdatfinl.functime()
                print 'Done in %.3g seconds.' % (timefinl - timeinit)
            
            if gdatfinl.verbtype > 0:
                print 'Accumulating lists...'
                timeinit = gdatfinl.functime()
            
            # lists of lists collected at each sample
            for strgvarb in gdatfinl.liststrgvarblistsamp:
                listtemp = [[[] for k in gdatfinl.indxproc] for j in gdatfinl.indxsamp]
                for j in gdatfinl.indxsamp:      
                    for k in gdatfinl.indxproc:
                        listtemp[j][k] = getattr(listgdatmodi[k], 'list' + strgpdfn + strgvarb)[j]
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listtemp)
            
            if gdatfinl.verbtype > 0:
                timefinl = gdatfinl.functime()
                print 'Done in %.3g seconds.' % (timefinl - timeinit)
            
            if not booltile:
                ## np.maximum likelihood sample 
                gdatfinl.maxmllikproc = np.empty(gdatfinl.numbproc)
                gdatfinl.indxswepmaxmllikproc = np.empty(gdatfinl.numbproc, dtype=int)
                gdatfinl.sampmaxmllikproc = np.empty((gdatfinl.numbproc, gdatfinl.fittnumbpara))
                for k in gdatfinl.indxproc:
                    gdatfinl.maxmllikproc[k] = listgdatmodi[k].maxmllikswep
                    gdatfinl.indxswepmaxmllikproc[k] = listgdatmodi[k].indxswepmaxmllik
                    gdatfinl.sampmaxmllikproc[k] = listgdatmodi[k].sampmaxmllik
            
                listsamp = getattr(gdatfinl, 'list' + strgpdfn + 'samp')
                listsampunit = getattr(gdatfinl, 'list' + strgpdfn + 'sampunit')

                # Gelman-Rubin test
                if gdatfinl.numbproc > 1:
                    if gdatfinl.verbtype > 0:
                        print 'Computing the Gelman-Rubin TS...'
                        timeinit = gdatfinl.functime()
                    gdatfinl.gmrbfixp = np.zeros(gdatfinl.fittnumbfixp)
                    gdatfinl.gmrbstat = np.zeros((gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbevtt))
                    for k in gdatfinl.fittindxfixp:
                        gdatfinl.gmrbfixp[k] = tdpy.mcmc.gmrb_test(listsamp[:, :, k])
                        if not np.isfinite(gdatfinl.gmrbfixp[k]):
                            gdatfinl.gmrbfixp[k] = 0.
                    listcntpmodl = getattr(gdatfinl, 'list' + strgpdfn + 'cntpmodl')
                    for i in gdatfinl.indxener:
                        for j in gdatfinl.indxpixl:
                            for m in gdatfinl.indxevtt:
                                gdatfinl.gmrbstat[i, j, m] = tdpy.mcmc.gmrb_test(listcntpmodl[:, :, i, j, m])
                    if gdatfinl.verbtype > 0:
                        timefinl = gdatfinl.functime()
                        print 'Done in %.3g seconds.' % (timefinl - timeinit)

                # calculate the autocorrelation of the chains
                if gdatfinl.verbtype > 0:
                    print 'Computing the autocorrelation of the chains...'
                    timeinit = gdatfinl.functime()
                gdatfinl.atcrcntp = np.empty((gdatfinl.numbproc, gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbevtt, gdatfinl.numbsamp / 2))
                gdatfinl.timeatcrcntp = np.empty((gdatfinl.numbproc, gdatfinl.numbener, gdatfinl.numbpixl, gdatfinl.numbevtt))
                gdatfinl.atcrpara = np.empty((gdatfinl.numbproc, gdatfinl.fittnumbpara, gdatfinl.numbsamp / 2))
                gdatfinl.timeatcrpara = np.empty((gdatfinl.numbproc, gdatfinl.fittnumbpara))
                for k in gdatfinl.indxproc:
                    gdatfinl.atcrpara[k, :, :], gdatfinl.timeatcrpara[k, :] = tdpy.mcmc.retr_timeatcr(listsamp[:, k, :], verbtype=gdatfinl.verbtype)
                    listcntpmodl = getattr(gdatfinl, 'list' + strgpdfn + 'cntpmodl')
                    gdatfinl.atcrcntp[k, :], gdatfinl.timeatcrcntp[k, :] = tdpy.mcmc.retr_timeatcr(listcntpmodl[:, k, :, :, :], verbtype=gdatfinl.verbtype)
                timeatcrcntpmaxm = np.amax(gdatfinl.timeatcrcntp)
                gdatfinl.timeatcrcntpmaxm = np.amax(timeatcrcntpmaxm)
                
                if gdatfinl.verbtype > 0:
                    timefinl = gdatfinl.functime()
                    print 'Done in %.3g seconds.' % (timefinl - timeinit)
                
                setattr(gdatfinl, 'list' + strgpdfn + 'sampproc', np.copy(getattr(gdatfinl, 'list' + strgpdfn + 'samp')))

            # flatten the list chains from different walkers
            for strgvarb in gdatfinl.liststrgvarblistsamp:
                listtemp = []
                listinpt = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                for j in gdatfinl.indxsamp:      
                    for k in gdatfinl.indxproc:
                        listtemp.append(listinpt[j][k])
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listtemp)
            
            # flatten the np.array chains from different walkers
            for strgvarb in gdatinit.liststrgvarbarryflat:
                inpt = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                shap = [inpt.shape[0] * inpt.shape[1]] + list(inpt.shape[2:])
                setattr(gdatfinl, 'list' + strgpdfn + strgvarb, inpt.reshape(shap))
            listsamp = getattr(gdatfinl, 'list' + strgpdfn + 'samp')
            listsampunit = getattr(gdatfinl, 'list' + strgpdfn + 'sampunit')
        
            if booltile:
                
                # temp
                #if n > 5:
                #    print 'Quitting...'
                #    break

                liststrgtile.append(rtagmodi.split('_')[-2][-4:])
                listrtaggood.append(rtagmodi)
                indxrtaggood.append(n)
                indxtiletemp += 1
                
                if len(liststrgtile) == 1:
                    for strgfeat in gdatfinl.refrliststrgfeattotl:
                        refrfeattile = [[] for q in gdatfinl.indxrefr]
                        setattr(gdatfinl, 'refr' + strgfeat, refrfeattile)
                
                    for strgvarb in gdatfinl.liststrgvarbarrysamp:
                        if not strgvarb in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                            listvarb = []
                            setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listvarb)
                        else:
                            hist = np.zeros_like(getattr(listgdatmodi[0], 'list' + strgpdfn + strgvarb))
                            setattr(gdatfinl, 'list' + strgpdfn + strgvarb, hist)
                
                    for name, valu in gdatfinl.__dict__.iteritems():
                        if name.startswith('refrhist'):
                            setattr(gdatfinl, name, np.zeros_like(getattr(gdatfinl, name)))
                            
                #for strgfeat in gdatfinl.refrliststrgfeattotl:
                #    refrfeattile = getattr(gdatfinl, 'refr' + strgfeat)
                #    #refrfeat = getattr(gdatfinl, 'refr' + strgfeat)
                #    print 'refrfeat'
                #    print refrfeat
                #    print
                #    refrfeat = [[] for q in gdatfinl.indxrefr]
                #    for q in gdatfinl.indxrefr:
                #        if strgfeat in gdatfinl.refrliststrgfeat[q]:
                #            refrfeat[q].append(refrfeattile[q])
                
                for strgvarb in gdatfinl.liststrgvarbarrysamp:
                    if strgvarb in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                        # temp
                        if 'spec' in strgvarb:
                            continue
                        hist = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                        hist += getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                    
                for name, valu in gdatfinl.__dict__.iteritems():
                    if name.startswith('refrhist'):
                        hist = getattr(gdatfinl, name)
                        hist += getattr(gdatfinl, name)

                print 'Done with the tile number %d, run number %d...' % (indxtiletemp, n)
                print
        
        if booltile:
            gdatfinl.pathplotrtag = gdatfinl.pathimag + rtagfinl + '/'
            make_fold(gdatfinl)
            indxrtaggood = np.array(indxrtaggood).astype(int)
            numbrtaggood = indxrtaggood.size
            numbtile = numbrtaggood
            print 'Found %d tiles with run tags:' % numbrtaggood
            for indxrtaggoodtemp in indxrtaggood:
                print rtag[indxrtaggoodtemp]
            print

            # np.concatenate reference elements from different tiles
            #for strgfeat in gdatfinl.refrliststrgfeattotl:
            #    refrfeat = getattr(gdatfinl, 'refr' + strgfeat, refrfeat)
            #    for q in gdatfinl.indxrefr:
            #        if strgfeat in gdatfinl.refrliststrgfeat[q]:
            #            refrfeat[q] = np.concatenate(refrfeat[q], axis=1)
            
            for strgvarb in gdatfinl.liststrgvarbarrysamp:
                if not strgvarb in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                    listvarb = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                    if 'assc' in strgvarb:
                        numbrefrelemtotl = 0
                        for k, varbrsmp in enumerate(listvarb):
                            numbrefrelemtotl += varbrsmp.shape[1]
                        shap = [gdatfinl.numbsamptotl, numbrefrelemtotl]
                        listvarbtemp = np.empty(shap)
                        cntr = 0
                        for k, varb in enumerate(listvarb):
                            listvarbtemp[:, cntr:cntr+varb.shape[1]] = varb
                            cntr += varb.shape[1]
                    else:
                        shap = [gdatfinl.numbsamptotl * numbtile] + list(listvarb[0].shape[1:])
                        listvarbtemp = np.empty(shap)
                        for k, varb in enumerate(listvarb):
                            listvarbtemp[k*gdatfinl.numbsamptotl:(k+1)*gdatfinl.numbsamptotl, ...] = varb
                    setattr(gdatfinl, 'list' + strgpdfn + strgvarb, listvarbtemp)
        else:
            # np.maximum likelihood sample
            if gdatfinl.fittnumbtrap > 0:
                listindxelemfull = getattr(gdatfinl, 'list' + strgpdfn + 'indxelemfull')
            listllik = getattr(gdatfinl, 'list' + strgpdfn + 'llik')
            listlliktotl = getattr(gdatfinl, 'list' + strgpdfn + 'lliktotl')
            indxsamptotlmlik = np.argmax(np.sum(np.sum(np.sum(listllik, 3), 2), 1))
            
            # copy the np.maximum likelihood sample
            for strgvarb in gdatfinl.liststrgvarbarrysamp:
                setattr(gdatfinl, 'mlik' + strgvarb, getattr(gdatfinl, 'list' + strgpdfn + strgvarb)[indxsamptotlmlik, ...])
            for strgvarb in gdatfinl.liststrgvarblistsamp:
                setattr(gdatfinl, 'mlik' + strgvarb, getattr(gdatfinl, 'list' + strgpdfn + strgvarb)[indxsamptotlmlik])

            # temp -- dont gdatfinl.listllik and gdatfinl.listsamp have the same dimensions?
            gdatfinl.mliksamp = getattr(gdatfinl, 'list' + strgpdfn + 'samp')[indxsamptotlmlik, :]
            gdatfinl.mliksamp = getattr(gdatfinl, 'list' + strgpdfn + 'samp')[indxsamptotlmlik, :]
            #if gdatfinl.fittnumbtrap > 0:
            #    gdatfinl.mlikindxelemfull = listindxelemfull[indxsamptotlmlik]
            gdatfinl.mlikfixp = gdatfinl.mliksamp[gdatfinl.fittindxfixp]
            for k, namefixp in enumerate(gdatfinl.fittnamefixp):
                setattr(gdatfinl, 'mlik' + namefixp, gdatfinl.mlikfixp[k])

            # add execution times to the chain output
            gdatfinl.timereal = np.zeros(gdatfinl.numbproc)
            gdatfinl.timeproc = np.zeros(gdatfinl.numbproc)
            for k in gdatfinl.indxproc:
                gdatfinl.timereal[k] = listgdatmodi[k].timereal
                gdatfinl.timeproc[k] = listgdatmodi[k].timeproc
        
            # find the np.maximum likelihood and posterior over the chains
            gdatfinl.indxprocmaxmllik = np.argmax(gdatfinl.maxmllikproc)
            #gdatfinl.maxmlliktotl = gdatfinl.maxmllikproc[gdatfinl.indxprocmaxmllik]
            gdatfinl.indxswepmaxmllik = gdatfinl.indxprocmaxmllik * gdatfinl.numbsamp + gdatfinl.indxswepmaxmllikproc[gdatfinl.indxprocmaxmllik]
            gdatfinl.sampmaxmllik = gdatfinl.sampmaxmllikproc[gdatfinl.indxprocmaxmllik, :]
                
            if strgpdfn == 'post':
                levipost = retr_levipost(listlliktotl)
                setattr(gdatfinl, strgpdfn + 'levipost', levipost)
            
            if strgpdfn == 'prio':
                leviprio = np.log(np.mean(np.exp(listlliktotl)))
                setattr(gdatfinl, strgpdfn + 'leviprio', leviprio)
            
        # parse the sample vector
        listfixp = listsamp[:, gdatfinl.fittindxfixp]
        for k, namefixp in enumerate(gdatfinl.fittnamefixp):
            setattr(gdatfinl, 'list' + strgpdfn + namefixp, listfixp[:, k])
        setattr(gdatfinl, 'list' + strgpdfn + 'fixp', listfixp)

        if strgpdfn == 'post' and gdatfinl.checprio:
            pathoutprtag = retr_pathoutprtag(rtag)
            path = pathoutprtag + 'gdatfinlprio'
            try:
                gdatprio = readfile(path)
            except:
                proc_finl(gdat=gdatfinl, strgpdfn='prio', listnamevarbproc=listnamevarbproc, forcplot=forcplot)
        else:
            gdatprio = None
        
        # post process samples
        ## bin element features
        if gdatfinl.verbtype > 0:
            print 'Binning the probabilistic catalog spatially...'
            timeinit = gdatfinl.functime()
        
        if not booltile:
            if gdatfinl.fittnumbtrap > 0:
                if gdatfinl.fittboolelemspatanyy:
                    histlgalbgalelemstkd = [[] for l in gdatfinl.fittindxpopl]
                
                    listlgal = getattr(gdatfinl, 'list' + strgpdfn + 'lgal')
                    listbgal = getattr(gdatfinl, 'list' + strgpdfn + 'bgal')
                    for l in gdatfinl.fittindxpopl:
                        if gdatfinl.fittelemtype[l] != 'lghtline':
                            numb = len(gdatfinl.fittliststrgfeatsign[l])
                            histlgalbgalelemstkd[l] = np.zeros((gdatfinl.numbbgalpntsprob, gdatfinl.numblgalpntsprob, gdatfinl.numbbinsplot, numb))
                            temparry = np.concatenate([listlgal[n][l] for n in gdatfinl.indxsamptotl])
                            temp = np.empty((len(temparry), 3))
                            temp[:, 0] = temparry
                            temp[:, 1] = np.concatenate([listbgal[n][l] for n in gdatfinl.indxsamptotl])
                            for k, strgfeat in enumerate(gdatfinl.fittliststrgfeatsign[l]):
                                temp[:, 2] = np.concatenate([getattr(gdatfinl, 'list' + strgpdfn + strgfeat)[n][l] for n in gdatfinl.indxsamptotl])
                                bins = getattr(gdatfinl, 'bins' + strgfeat)
                                histlgalbgalelemstkd[l][:, :, :, k] = histogramdd(temp, bins=(gdatfinl.binslgalpntsprob, gdatfinl.binsbgalpntsprob, bins))[0]
                    setattr(gdatfinl, strgpdfn + 'histlgalbgalelemstkd', histlgalbgalelemstkd)

            if gdatfinl.verbtype > 0:
                timefinl = gdatfinl.functime()
                print 'Done in %.3g seconds.' % (timefinl - timeinit)

            ## construct a condensed catalog of elements
            if gdatfinl.boolcondcatl and gdatfinl.fittnumbtrap > 0:
                
                if gdatfinl.verbtype > 0:
                    print 'Constructing a condensed catalog...'
                    timeinit = gdatfinl.functime()
                
                retr_condcatl(gdat)
            
                if gdatfinl.verbtype > 0:
                    timefinl = gdatfinl.functime()
                    print 'Done in %.3g seconds.' % (timefinl - timeinit)

            # construct lists of samples for each proposal type
            listindxproptype = getattr(gdatfinl, 'list' + strgpdfn + 'indxproptype')
            listboolpropaccp = getattr(gdatfinl, 'list' + strgpdfn + 'boolpropaccp')
            listboolpropfilt = getattr(gdatfinl, 'list' + strgpdfn + 'boolpropfilt')
            listindxsamptotlproptotl = []
            listindxsamptotlpropfilt = []
            listindxsamptotlpropaccp = []
            listindxsamptotlpropreje = []
            for n in gdatfinl.indxproptype:
                listindxsamptotlproptotl.append(np.where(listindxproptype == gdatfinl.indxproptype[n])[0])
                listindxsamptotlpropaccp.append(np.intersect1d(np.where(listindxproptype == gdatfinl.indxproptype[n])[0], np.where(listboolpropaccp)[0]))
                listindxsamptotlpropfilt.append(np.intersect1d(np.where(listindxproptype == gdatfinl.indxproptype[n])[0], np.where(listboolpropfilt)[0]))
                listindxsamptotlpropreje.append(np.intersect1d(np.where(listindxproptype == gdatfinl.indxproptype[n])[0], np.where(np.logical_not(listboolpropaccp))[0]))
                if listindxsamptotlproptotl[n].size == 0:
                    accp = 0.
                else:
                    accp = float(listindxsamptotlpropaccp[n].size) / listindxsamptotlproptotl[n].size
                setattr(gdatfinl, 'accp' + gdatfinl.nameproptype[n], accp)

            setattr(gdatfinl, 'list' + strgpdfn + 'indxsamptotlproptotl', listindxsamptotlproptotl)
            setattr(gdatfinl, 'list' + strgpdfn + 'indxsamptotlpropaccp', listindxsamptotlpropaccp)
            setattr(gdatfinl, 'list' + strgpdfn + 'indxsamptotlpropreje', listindxsamptotlpropreje)
       
        if gdatfinl.fittnumbtrap > 0 and strgpdfn == 'post':
            if gdatfinl.datatype == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:
                    if gdatfinl.rtagmock is not None:
                        path = gdatfinl.pathoutprtagmock + 'gdatfinlpost'
                        gdatmock = readfile(path)
                    
        # posterior corrections
        if gdatfinl.fittnumbtrap > 0 and strgpdfn == 'post':

            ## perform corrections
            if gdatfinl.datatype == 'inpt':
                if gdatfinl.boolcrex or gdatfinl.boolcrin:

                    for liststrgcompvarbhist in gdatfinl.liststrgvarbhist:
                        strgvarb = liststrgcompvarbhist[0]

                        if liststrgcompvarbhist[1].startswith('aerr') or len(liststrgcompvarbhist[2]) > 0 and liststrgcompvarbhist[2].startswith('aerr'):
                            continue
                        if liststrgcompvarbhist[1] == 'spec' or liststrgcompvarbhist[1] == 'deflprof' or liststrgcompvarbhist[1] == 'specplot':
                            continue
                        if len(liststrgcompvarbhist[2]) > 0 and (liststrgcompvarbhist[2] == 'spec' or \
                                    liststrgcompvarbhist[2] == 'deflprof' or liststrgcompvarbhist[2] == 'specplot'):
                            continue
                        
                        ## internal correction
                        listhist = getattr(gdatfinl, 'list' + strgpdfn + strgvarb)
                        
                        for qq in gdatmock.indxrefr:
                            print 'liststrgcompvarbhist[3][qq]'
                            print liststrgcompvarbhist[3][qq]
                            l = int(liststrgcompvarbhist[3][qq].split('pop')[1][0])
                            qq = int(liststrgcompvarbhist[3][qq].split('pop')[2][0])
                            if liststrgcompvarbhist[1][-4:] in gdatfinl.listnamerefr and \
                                    (len(liststrgcompvarbhist[2]) == 0 or liststrgcompvarbhist[2][-4:] in gdatfinl.listnamerefr):
                                listhistincr = listhist
                            else:
                                if liststrgcompvarbhist[1][-4:] in gdatfinl.listnamerefr and len(liststrgcompvarbhist[2]) > 0:
                                    listcmpltrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostcmpl' + liststrgcompvarbhist[2] + 'pop%dpop%d' % (l, qq))], 2)
                                    listfdistrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostfdis' + liststrgcompvarbhist[2] + 'pop%dpop%d' % (qq, l))], 2)
                                elif len(liststrgcompvarbhist[2][:-4]) > 0 and liststrgcompvarbhist[2][-4:] in gdatfinl.listnamerefr:
                                    listcmpltrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostcmpl' + liststrgcompvarbhist[1] + 'pop%dpop%d' % (l, qq))], 1)
                                    listfdistrue = np.stack(gdatfinl.numbbinsplot * \
                                                    [getattr(gdatmock, 'listpostfdis' + liststrgcompvarbhist[1] + 'pop%dpop%d' % (qq, l))], 1)
                                else:
                                    listcmpltrue = getattr(gdatmock, 'listpostcmpl' + liststrgcompvarbhist[3][qq])
                                    listfdistrue = getattr(gdatmock, 'listpostfdis' + liststrgcompvarbhist[3][qq])
                                if len(liststrgcompvarbhist[2]) == 0:
                                    listcmplboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot))
                                    listfdisboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot))
                                    listhistboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot))
                                    for k in gdatfinl.indxbinsplot:
                                        listcmplboot[:, k] = choice(listcmpltrue[:, k], size=gdatfinl.numbsampboot)
                                        listfdisboot[:, k] = choice(listfdistrue[:, k], size=gdatfinl.numbsampboot)
                                        listhistboot[:, k] = choice(listhist[:, k], size=gdatfinl.numbsampboot)
                                else:
                                    listcmplboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot, gdatfinl.numbbinsplot))
                                    listfdisboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot, gdatfinl.numbbinsplot))
                                    listhistboot = np.empty((gdatfinl.numbsampboot, gdatfinl.numbbinsplot, gdatfinl.numbbinsplot))
                                    for a in gdatfinl.indxbinsplot:
                                        for b in gdatfinl.indxbinsplot:
                                            listcmplboot[:, a, b] = choice(listcmpltrue[:, a, b], size=gdatfinl.numbsampboot)
                                            listfdisboot[:, a, b] = choice(listfdistrue[:, a, b], size=gdatfinl.numbsampboot)
                                            listhistboot[:, a, b] = choice(listhist[:, a, b], size=gdatfinl.numbsampboot)
                                indxbadd = np.where(listcmplboot == -1)
                                indxbaddzero = np.where(listcmplboot == 0.)
                                listhistincr = listhistboot / listcmplboot * (1. - listfdisboot)
                                listhistincr[indxbadd] = -1.5
                                listhistincr[indxbaddzero] = 1.5
                            
                            gdatfinl.liststrgchan += ['incr' + liststrgcompvarbhist[4][qq]]
                            setattr(gdatfinl, 'listpostincr' + liststrgcompvarbhist[4][qq], listhistincr)
                        
                            ## external correction
                            for q in gdatfinl.indxrefr:
                                nametemp = liststrgcompvarbhist[1] 
                                if len(liststrgcompvarbhist[2]) > 0:
                                    nametemp += liststrgcompvarbhist[2]
                                nametemp += 'pop%dpop%dpop%d' % (q, qq, l)
                                crexhist = getattr(gdatfinl, 'crex' + nametemp)
                                if crexhist is not None:
                                    
                                    listhistexcr = listhistincr * crexhist 
                                    
                                    if crexhist.ndim == 1 and listhistincr.ndim == 3:
                                        raise Exception('')
                                    
                                    gdatfinl.liststrgchan += ['excr' + nametemp]
                                    setattr(gdatfinl, 'listpostexcr' + nametemp, listhistexcr)
                            
        # compute credible intervals
        if gdatfinl.verbtype > 0:
            print 'Computing credible intervals...'
            timeinit = gdatfinl.functime()
       
        for strgchan in gdatfinl.liststrgchan:
            
            if booltile:
                if strgchan in gdatfinl.liststrgvarbarryswep or strgchan in gdatfinl.liststrgvarblistsamp:
                    continue
                if not (strgchan.startswith('hist') or strgchan.startswith('incr') or strgchan.startswith('excr')):
                    continue

            if gdatfinl.fittnumbtrap > 0 and strgchan in [strgvarbhist[0] for strgvarbhist in gdatfinl.liststrgvarbhist]:
                if 'spec' in strgchan:
                    continue
            if strgchan == 'spec':
                continue

            listtemp = getattr(gdatfinl, 'list' + strgpdfn + strgchan)
            
            if isinstance(listtemp, list):
            
                if booltile:
                    continue

                # ensure that transdimensional lists are not included
                # temp
                if strgchan in gdatfinl.fittliststrgfeattotl or strgchan == 'indxelemfull':
                    continue

                pctltemp = []
                pmeatemp = []
                meditemp = []
                errrtemp = []
                stdvtemp = []
                numb = len(listtemp[0])
                
                for k in range(numb):
                    if isinstance(listtemp[0][k], list):
                        continue
                    shap = [gdatfinl.numbsamptotl] + list(listtemp[0][k].shape)
                    temp = np.zeros(shap)
                    for n in gdatfinl.indxsamptotl:
                        temp[n, ...] = listtemp[n][k]
                    
                    pctltempsing = tdpy.util.retr_pctlvarb(temp)
                    pmeatempsing = np.mean(temp, axis=0)
                    meditempsing = pctltempsing[0, ...]
                    errrtempsing = tdpy.util.retr_errrvarb(pctltempsing)
                    stdvtempsing = np.std(temp)
                    
                    pctltemp.append(pctltempsing)
                    pmeatemp.append(pmeatempsing)
                    meditemp.append(meditempsing)
                    errrtemp.append(errrtempsing)
                    stdvtemp.append(stdvtempsing)
            else:
                # this is needed for finding posterior moments of features of associated reference elements
                if 'asscref' in strgchan:
                    if listtemp.ndim != 2:
                        raise Exception('')
                    pmeatemp = np.zeros(listtemp.shape[1])
                    pctltemp = np.zeros([3] + [listtemp.shape[1]])
                    # temp -- this only works for 2D listtemp
                    for k in range(listtemp.shape[1]):
                        indxassc = np.where(np.isfinite(listtemp[:, k]))[0]
                        if indxassc.size > 0:
                            pctltemp[:, k] = tdpy.util.retr_pctlvarb(listtemp[indxassc, k])
                            pmeatemp[k] = np.mean(listtemp[indxassc, k])
                else:
                    pctltemp = tdpy.util.retr_pctlvarb(listtemp)
                    pmeatemp = np.mean(listtemp, axis=0)
                
                errrtemp = tdpy.util.retr_errrvarb(pctltemp)
                stdvtemp = np.std(pctltemp, axis=0)
                meditemp = pctltemp[0, ...]
                
                if strgchan in gdatfinl.listnamevarbcpct:
                    cpcttemp = np.empty([gdatfinl.numbsampcpct] + [3] + list(listtemp.shape[1:]))
                    for n in gdatfinl.indxsampcpct:
                        cpcttemp[n, ...] = tdpy.util.retr_pctlvarb(listtemp[:n+1, ...])
            
            setattr(gdatfinl, 'pctl' + strgpdfn + strgchan, pctltemp)
            setattr(gdatfinl, 'medi' + strgpdfn + strgchan, meditemp)
            setattr(gdatfinl, 'pmea' + strgpdfn + strgchan, pmeatemp)
            setattr(gdatfinl, 'errr' + strgpdfn + strgchan, errrtemp)
            setattr(gdatfinl, 'stdv' + strgpdfn + strgchan, stdvtemp)
            if strgchan in gdatfinl.listnamevarbcpct:
                setattr(gdatfinl, 'cpct' + strgpdfn + strgchan, cpcttemp)
        
        if not booltile:
            pmealliktotl = getattr(gdatfinl, 'pmea' + strgpdfn + 'lliktotl')
            stdvlliktotl = getattr(gdatfinl, 'stdv' + strgpdfn + 'lliktotl')
            minmlliktotl = np.amin(listlliktotl)
            maxmlliktotl = np.amax(listlliktotl)
            skewlliktotl = np.mean(((listlliktotl - pmealliktotl) / stdvlliktotl)**3)
            kurtlliktotl = np.mean(((listlliktotl - pmealliktotl) / stdvlliktotl)**4)
            setattr(gdatfinl, 'minm' + strgpdfn + 'lliktotl', minmlliktotl)
            setattr(gdatfinl, 'maxm' + strgpdfn + 'lliktotl', maxmlliktotl)
            setattr(gdatfinl, 'skew' + strgpdfn + 'lliktotl', skewlliktotl)
            setattr(gdatfinl, 'kurt' + strgpdfn + 'lliktotl', kurtlliktotl)

            if strgpdfn == 'post':
                infopost = retr_infofromlevi(pmealliktotl, levipost)
                setattr(gdatfinl, strgpdfn + 'infopost', infopost)
            if strgpdfn == 'post' and gdatfinl.checprio:
                leviprio = getattr(gdatprio, 'prioleviprio')
                infoprio = retr_infofromlevi(pmealliktotl, leviprio)
                setattr(gdatfinl, strgpdfn + 'infoprio', infoprio)
            
            bcom = maxmlliktotl - pmealliktotl
            setattr(gdatfinl, strgpdfn + 'bcom', bcom)
        
        listnametemp = ['lliktotl']
        if gdat.fittnumbtrap > 0:
            listnametemp += ['lpripena']

        for namevarbscal in listnametemp:
            listtemp = getattr(gdatfinl, 'list' + strgpdfn + namevarbscal)
            minm = np.amin(listtemp)
            maxm = np.amax(listtemp)
            setattr(gdatfinl, 'minm' + namevarbscal, minm)
            setattr(gdatfinl, 'maxm' + namevarbscal, maxm)
            setattr(gdatfinl, 'scal' + namevarbscal, 'self')
            retr_axis_wrap(gdatfinl, namevarbscal)
        
        if gdatfinl.checprio:
            for strgvarb in gdatfinl.listnamevarbscal:
                setp_pdfnvarb(gdatfinl, strgpdfn, strgvarb, strgvarb)
            for l0 in gdatfinl.fittindxpopl:
                for strgfeatfrst in gdatfinl.fittliststrgfeat[l0]:
                    if strgfeatfrst == 'spec' or strgfeatfrst == 'deflprof' or strgfeatfrst == 'specplot':
                        continue
                    setp_pdfnvarb(gdatfinl, strgpdfn, strgfeatfrst, 'hist' + strgfeatfrst + 'pop%d' % l0)
                    for strgfeatseco in gdatfinl.fittliststrgfeat[l0]:
                        if strgfeatseco == 'spec' or strgfeatseco == 'deflprof' or strgfeatseco == 'specplot':
                            continue
                        
                        if not checstrgfeat(strgfeatfrst, strgfeatseco):
                            continue
                                    
                        setp_pdfnvarb(gdatfinl, strgpdfn, strgfeatfrst, 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l0, nameseco=strgfeatseco)

            # calculate information gain
            if strgpdfn == 'post':
                for namevarbscal in gdatfinl.listnamevarbscal:
                    setp_info(gdatfinl, gdatprio, namevarbscal, namevarbscal)
                for l0 in gdatfinl.fittindxpopl:
                    for strgfeatfrst in gdatfinl.fittliststrgfeat[l0]:
                        if strgfeatfrst == 'spec' or strgfeatfrst == 'deflprof' or strgfeatfrst == 'specplot':
                            continue
                        setp_info(gdatfinl, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + 'pop%d' % l0)
                        for strgfeatseco in gdatfinl.fittliststrgfeat[l0]:
                            if strgfeatseco == 'spec' or strgfeatseco == 'deflprof' or strgfeatseco == 'specplot':
                                continue
                            
                            if not checstrgfeat(strgfeatfrst, strgfeatseco):
                                continue
                                    
                            setp_info(gdatfinl, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l0, nameseco=strgfeatseco)

        if gdatfinl.verbtype > 0:
            timefinl = gdatfinl.functime()
            print 'Done in %.3g seconds.' % (timefinl - timeinit)
        
        # flatten the np.arrays which have been collected at each sweep
        #setattr(gdat, 'list' + strgpdfn + strgpdfntemp + 'flat', getattr(gdat, 'list' + strgpdfn + strgpdfntemp + 'totl').flatten())
        if not booltile:
            # memory usage
            listmemoresi = getattr(gdatfinl, 'list' + strgpdfn + 'memoresi')
            gdatfinl.meanmemoresi = np.mean(listmemoresi, 1)
            gdatfinl.derimemoresi = (gdatfinl.meanmemoresi[-1] - gdatfinl.meanmemoresi[0]) / gdatfinl.numbswep

            gdatfinl.timerealtotl = time.time() - gdatfinl.timerealtotl
            gdatfinl.timeproctotl = time.clock() - gdatfinl.timeproctotl
            gdatfinl.timeproctotlswep = gdatfinl.timeproctotl / gdatfinl.numbswep
            
            if gdatfinl.timeatcrcntpmaxm == 0.:
                gdatfinl.timeprocnorm = 0.
            else:
                gdatfinl.timeprocnorm = gdatfinl.timeproctotlswep / gdatfinl.timeatcrcntpmaxm
   
        # write the final gdat object
        path = gdatfinl.pathoutprtag + 'gdatfinl' + strgpdfn

        if gdatfinl.verbtype > 0:
            print 'Writing gdatfinl to %s...' % path
        writfile(gdatfinl, path) 
       
        filestat = open(gdatfinl.pathoutprtag + 'stat.txt', 'a')
        filestat.write('gdatfinl%s written.\n' % strgpdfn)
        filestat.close()
   
        if not booltile:
            if gdatfinl.verbtype > 0:
                for k in gdatfinl.indxproc:
                    print 'Process %d has been completed in %d real seconds, %d CPU seconds.' % (k, gdatfinl.timereal[k], gdatfinl.timeproc[k])
                print 'Parent process has run in %d real seconds, %d CPU seconds.' % (gdatfinl.timerealtotl, gdatfinl.timeproctotl)

    print 'checking plotfinl%s...' % strgpdfn


    booltemp = chec_statfile(rtagfinl, 'plotfinl', strgpdfn)
    if booltemp:
        print 'Final plots already written.'
    else:
        if strgpdfn == 'post' and gdatfinl.checprio:
            path = pathoutprtag + 'gdatfinlprio'
            gdatprio = readfile(path)
        else:
            gdatprio = None
        
        if gdatfinl.makeplot and getattr(gdatfinl, 'makeplotfinl' + strgpdfn) or forcplot:
            plot_finl(gdatfinl, gdatprio=gdatprio, strgpdfn=strgpdfn, gdatmock=gdatmock, booltile=booltile)
            filestat = open(gdatfinl.pathoutprtag + 'stat.txt', 'a')
            filestat.write('plotfinl%s written.\n' % strgpdfn)
            filestat.close()
    print


def retr_axis_wrap(gdat, namevarbscal):

    minm = getattr(gdat, 'minm' + namevarbscal)
    maxm = getattr(gdat, 'maxm' + namevarbscal)
    if namevarbscal in gdat.fittnamefixp:
        scal = getattr(gdat, 'fittscal' + namevarbscal)
    else:
        scal = getattr(gdat, 'scal' + namevarbscal)
    if gdat.booldiagmode:
        if not isinstance(scal, str):
            print 'namevarbscal'
            print namevarbscal
            print 'scal'
            print scal
            raise Exception('Parameter scaling is bad.')
    retr_axis(gdat, namevarbscal, minm, maxm, gdat.numbbinspdfn, scal=scal)


def retr_listgdat(listrtag, typegdat='finlpost'):
   
    listgdat = []
    for rtag in listrtag:
        pathoutprtag = retr_pathoutprtag(rtag)
        path = pathoutprtag + 'gdat%s' % typegdat
        listgdat.append(readfile(path))

    return listgdat


def make_fold(gdat):

    for strgpdfn in gdat.liststrgpdfn:
        setattr(gdat, 'path' + strgpdfn, gdat.pathplotrtag + strgpdfn + '/') 
        path = getattr(gdat, 'path' + strgpdfn)

        for nameseco in ['finl', 'fram', 'anim', 'opti']:
            setattr(gdat, 'path' + strgpdfn + nameseco, path + nameseco + '/')
        
        for nameseco in ['diag', 'lpac', 'varbscal', 'cond', 'varbscalproc']:
            setattr(gdat, 'path' + strgpdfn + 'finl' + nameseco, path + 'finl/' + nameseco + '/')
        
        for n in gdat.indxproptype:
            setattr(gdat, 'path' + strgpdfn + 'finl' + gdat.nameproptype[n], path + 'finl/lpac/' + gdat.nameproptype[n] + '/')

        for namethrd in ['hist', 'trac', 'join', 'cova']:
            setattr(gdat, 'path' + strgpdfn + 'finlvarbscal' + namethrd, path + 'finl/varbscal/' + namethrd + '/')
            
        for strgphas in gdat.liststrgphas + ['init']:
            liststrgfold = getattr(gdat, 'liststrgfold' + strgphas)
            for nameseco in liststrgfold:
                if strgphas == 'init':
                    if nameseco == 'assc/' or nameseco.startswith('cmpl') or nameseco.startswith('fdis'):
                        continue
                    setattr(gdat, 'path' + strgphas + nameseco[:-1], gdat.pathplotrtag + 'init/' + nameseco + '/')
                else:
                    setattr(gdat, 'path' + strgpdfn + strgphas + nameseco[:-1], path + strgphas + '/' + nameseco + '/')
    gdat.pathinfo = gdat.pathplotrtag + 'info/'
    
    ## make the directories 
    for attr, valu in gdat.__dict__.iteritems():
        if attr.startswith('path'):
            os.system('mkdir -p %s' % valu)


def make_cmapdivg(strgcolrloww, strgcolrhigh):
    
    funccolr = mpl.colors.ColorConverter().to_rgb
    
    colrloww = funccolr(strgcolrloww)
    colrhigh = funccolr(strgcolrhigh)
    
    cmap = make_cmap([colrloww, funccolr('white'), 0.5, funccolr('white'), colrhigh])

    return cmap


def make_cmap(seq):
    
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    
    return mpl.colors.LinearSegmentedColormap('CustomMap', cdict)


def setp_pdfnvarb(gdat, strgpdfn, name, namefull, nameseco=None):
    
    meanvarb = getattr(gdat, 'mean' + name)
    listvarb = getattr(gdat, 'list' + strgpdfn + namefull)
    
    if listvarb.ndim == 1:
        shaptemp = [gdat.numbbinspdfn, 1]
    else:
        shaptemp = [gdat.numbbinspdfn] + list(listvarb.shape[1:])
    pdfn = np.empty(shaptemp)
    if listvarb.ndim == 1:
        binsvarb = getattr(gdat, 'bins' + name)
        deltvarb = getattr(gdat, 'delt' + name)
        pdfn[:, 0] = np.histogram(listvarb, bins=binsvarb)[0].astype(float)
        pdfn[:, 0] /= np.sum(pdfn[:, 0])
        pdfn[:, 0] /= deltvarb
    else:
        binsvarb = np.linspace(0, gdat.fittmaxmnumbelemtotl, 51)
        
    if listvarb.ndim == 2:
        for k in range(listvarb.shape[1]):
            pdfn[:, k] = np.histogram(listvarb[:, k], bins=binsvarb)[0].astype(float)
            pdfn[:, k] /= np.sum(pdfn[:, k])
        pdfn *= 50.
    if listvarb.ndim == 3:
        for k in range(listvarb.shape[1]):
            for m in range(listvarb.shape[2]):
                pdfn[:, k, m] = np.histogram(listvarb[:, k, m], bins=binsvarb)[0].astype(float)
                pdfn[:, k, m] /= np.sum(pdfn[:, k, m])
        pdfn *= 2500.
    pdfn[np.where(pdfn < 1e-50)[0]] = 1e-50
    
    setattr(gdat, 'pdfn' + strgpdfn + namefull, pdfn)


def setp_info(gdat, gdatprio, name, namefull, nameseco=None, namesecofull=None):
    
    listpost = getattr(gdat, 'listpost' + namefull)
    listprio = getattr(gdatprio, 'listprio' + namefull)
    pdfnpost = getattr(gdat, 'pdfnpost' + namefull)
    pdfnprio = getattr(gdatprio, 'pdfnprio' + namefull)
    if listpost.ndim == 3:
        infodens = np.empty((gdat.numbbinspdfn, listpost.shape[1], listpost.shape[2]))
        info = np.empty((listpost.shape[1], listpost.shape[2]))
        pvks = np.empty((listpost.shape[1], listpost.shape[2]))
    else:
        if listpost.ndim == 1:
            numbtemp = 1
        else:
            numbtemp = listpost.shape[1]
        infodens = np.empty((gdat.numbbinspdfn, numbtemp))
        info = np.empty(numbtemp)
        pvks = np.empty(numbtemp)
    if listpost.ndim == 1:
        listpost = listpost[:, None]
        listprio = listprio[:, None]
        deltvarb = getattr(gdat, 'delt' + name)
    else:
        if listpost.ndim == 2:
            deltvarb = 1. / 50
        else:
            deltvarb = 1. / 50**2
    
    if listpost.ndim == 1 or listpost.ndim == 2:
        for k in range(listpost.shape[1]):
            infodens[:, k] = retr_infodens(pdfnpost[:, k], pdfnprio[:, k])
            info[k] = np.sum(infodens[:, k] * deltvarb)
            temp, pvks[k] = sp.stats.ks_2samp(listpost[:, k], listprio[:, k])
    if listpost.ndim == 3:
        for k in range(listpost.shape[1]):
            for m in range(listpost.shape[2]):
                infodens[:, k, m] = retr_infodens(pdfnpost[:, k, m], pdfnprio[:, k, m])
                info[k, m] = np.sum(infodens[:, k, m] * deltvarb)
                temp, pvks[k, m] = sp.stats.ks_2samp(listpost[:, k, m], listprio[:, k, m])
    
    setattr(gdat, 'pvks' + namefull, pvks)
    setattr(gdat, 'infodens' + namefull, infodens)
    setattr(gdat, 'info' + namefull, info)


def chec_statfile(rtag, strggdat, strgpdfn, verbtype=1):
    
    pathoutprtag = retr_pathoutprtag(rtag)
    
    # check the status file
    if not os.path.isfile(pathoutprtag + 'stat.txt'):
        if verbtype > 0:
            print 'pathoutprtag'
            print pathoutprtag
            print 'stat.txt not found.'
            print
        return False

    # check the global object
    filestat = open(pathoutprtag + 'stat.txt', 'r')
    booltemp = False
    for line in filestat:
        if line == strggdat + strgpdfn + ' written.\n':
            booltemp = True
    
    filestat.close()
    if not booltemp:
        if verbtype > 0:
            print 'bad %s status.' % (strggdat + strgpdfn)
            print
        return False
    else:
        return True


def retr_los3(dlos, lgal, bgal):

    dglc = np.sqrt(8.5e3**2 + dlos**2 - 2. * dlos * 8.5e3 * np.cos(bgal) * np.cos(lgal))
    thet = np.arccos(np.sin(bgal) * dlos / dglc)
    phii = np.arcsin(np.sqrt(np.cos(bgal)**2 * dlos**2 + 8.5e3**2 - 2 * dlos * np.cos(bgal) * 8.5e3) / dglc)
    
    return dglc, thet, phii


def retr_glc3(dglc, thet, phii):

    xpos = dglc * np.sin(thet) * np.cos(phii)
    ypos = dglc * np.sin(thet) * np.sin(phii)
    zpos = dglc * np.cos(thet)
    dlos = np.sqrt(zpos**2 + xpos**2 + (8.5e3 - ypos)**2)
    lgal = np.arctan2(8.5e3 - ypos, xpos) - np.pi / 2
    bgal = np.arcsin(zpos / dlos)
   
    return dlos, lgal, bgal


def retr_lumipuls(geff, magf, per0):

    # temp -- this is bolometric luminosity np.whereas dictelem[l]['flux'] is differential!
    lumi = 9.6e33 * (geff / 0.2) * (magf / 10**8.5)**2 * (3e-3 / per0)*4

    return lumi


def retr_lumi(gdat, flux, dlos, reds=None):

    lumi = flux * 4. * np.pi * dlos**2 * gdat.prsccmtr**2 / gdat.ergsgevv
    
    # temp
    # redshift correction
    if reds is not None:
        lumi *= (1. + reds)**2

    return lumi


def retr_flux(gdat, lumi, dlos, reds=None):

    flux = lumi / 4. / np.pi / dlos**2 / gdat.prsccmtr**2 * gdat.ergsgevv
    
    # temp
    # redshift correction
    if reds is not None:
        pass

    return flux


def retr_per1(per0, magf):

    per1 = 3.3e-20 * (magf / 10**8.5)**2 * (3e-3 / per0)

    return per1


def retr_dlosgalx(lgal, bgal, dglc):

    # temp -- this is obviously wrong
    dlos = 8.5e3 - dglc

    return dlos


def retr_arryfromlist(listtemp):
    
    shap = [len(listtemp)] + list(listtemp[0].shape)
    arry = np.empty(shap)
    for k in range(len(listtemp)):
        arry[k, ...] = listtemp[k]
    
    return arry


def proc_cntpdata(gdat):

    # exclude voxels with vanishing exposure
    ## data counts
    if gdat.datatype == 'inpt':
        gdat.cntpdata = retr_cntp(gdat, gdat.sbrtdata)
    
    # data variance
    gdat.varidata = np.maximum(gdat.cntpdata, 1.)

    # correct the likelihoods for the constant data dependent factorial
    gdat.llikoffs = -sp.special.gammaln(gdat.cntpdata + 1)

    ## spatial average
    gdat.sbrtdatamean, gdat.sbrtdatastdv = retr_spatmean(gdat, gdat.cntpdata, boolcntp=True)
    
    retr_datatick(gdat)
    
    # 1-point function of the data counts
    for m in gdat.indxevtt:
        if gdat.numbpixl > 1:
            for i in gdat.indxener: 
                histcntp = np.histogram(gdat.cntpdata[i, :, m], bins=gdat.binscntpdata)[0]
                setattr(gdat, 'histcntpdataen%02devt%d' % (i, m), histcntp)
        else:
            histcntp = np.histogram(gdat.cntpdata[:, 0, m], bins=gdat.binscntpdata)[0]
            setattr(gdat, 'histcntpdataevt%d' % m, histcntp)

    # obtain cartesian versions of the maps
    if gdat.pixltype == 'cart':
        ## data counts
        gdat.cntpdatacart = np.zeros((gdat.numbener, gdat.numbpixlcart, gdat.numbevtt))
        gdat.cntpdatacart[:, gdat.indxpixlrofi, :] = gdat.cntpdata
        gdat.cntpdatacart = gdat.cntpdatacart.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
   

def retr_infodens(pdfnpost, pdfnprio):
    
    infodens = pdfnpost * np.log(pdfnpost / pdfnprio)

    return infodens


def retr_llik(gdat, strgmodl, cntpmodl):
   
    if gdat.liketype == 'pois':
        llik = gdat.cntpdata * np.log(cntpmodl) - cntpmodl
    if gdat.liketype == 'gaus':
        llik = -0.5 * (gdat.cntpdata - cntpmodl)**2 / gdat.varidata
     
    return llik


def retr_mapsgaus(gdat, lgal, bgal, spec, size, ellp, angl):
    
    rttrmatr = np.array([[np.cos(angl), -np.sin(angl)], [np.sin(angl), np.cos(angl)]])
    icovmatr = np.array([[1. / ((1. - ellp) * size)**2, 0.], [0., 1. / size**2]])

    posi = np.array([lgalgrid - lgal, bgalgrid - bgal])
    mapsgaus = flux * np.exp(-0.5 * np.sum(posi * tensordot(self.icovmatr, posi, (1,0)), 0)) / size**2 / (1. - ellp)
        
    return mapsgaus


def retr_sbrtsers(gdat, lgalgrid, bgalgrid, lgal, bgal, spec, size, ellp, angl, seri=np.array([4.])):
   
    lgalrttr = (1. - ellp) * (np.cos(angl) * (lgalgrid - lgal) - np.sin(angl) * (bgalgrid - bgal))
    bgalrttr = np.sin(angl) * (lgalgrid - lgal) + np.cos(angl) * (bgalgrid - bgal) 
    angl = np.sqrt(lgalrttr**2 + bgalrttr**2)
    
    # interpolate pixel-convolved Sersic surface brightness
    if gdat.serstype == 'intp':

        shapinpt = angl.shape 
        inpt = np.empty(list(shapinpt) + [3])
        inpt[..., 0] = angl
        inpt[..., 1] = size
        inpt[..., 2] = seri
        
        sbrtsers = spec[:, None, None] * sp.interpolate.interpn((gdat.binslgalsers, gdat.binshalfsers, gdat.binsindxsers), gdat.sersprof, inpt)[None, :, None]
    
    # evaluate directly de Vaucouleurs
    if gdat.serstype == 'vauc':
        sbrtsers = spec[:, None, None] * retr_sbrtsersnorm(angl, size)[None, :, None]
    
    return sbrtsers


def retr_sbrtsersnorm(angl, halfsers, indxsers=4.):

    ## this approximation works for 0.5  < indx < 10
    factsers = 1.9992 * indxsers - 0.3271
    
    ## surface brightness profile at the half-light radius for a 1 erg cm^-2 s^-1 A^-1 source
    if indxsers == 4.:
        sbrthalf = 1. / 7.2 / np.pi / halfsers**2
    else:
        sbrthalf = 1. / 2. / np.pi / np.exp(factsers) * factsers**(2 * indxsers) / indxsers / sp.special.gamma(2. * indxsers) / halfsers**2
                
    ## surface brightness profile
    sbrtsers = sbrthalf * np.exp(-factsers * ((angl / halfsers)**(1. / indxsers) - 1.))
    
    return sbrtsers


def copytdgu(varb):
    
    if isinstance(varb, np.ndarray):
        return np.copy(varb)
    else:
        return deepcopy(varb)


def prep_gdatmodi(gdat, gdatmodi, gdatobjt, strgmodl):
    
    if gdat.verbtype > 1:
        print 'prep_gdatmodi()'
    
    strgpfixthis = retr_strgpfix('this', strgmodl)
    for k, name in enumerate(gdat.listnamechro):
        setattr(gdatmodi, 'thischro' + name, 0.)
    for namevarbstat in gdat.listnamevarbstat:
        temp = copytdgu(getattr(gdatobjt, strgpfixthis + namevarbstat))
        setattr(gdatmodi, strgpfixthis + namevarbstat, temp)
    

def proc_anim(rtag):
    
    pathoutprtag = retr_pathoutprtag(rtag)
    
    print 'Making animations of frame plots for %s...' % rtag
    
    path = pathoutprtag + 'gdatinit'
    gdat = readfile(path)
    for strgpdfn in gdat.liststrgpdfn:
        for nameextn in gdat.liststrgfoldanim:
            
            pathframextn = gdat.pathimag + rtag + '/' + strgpdfn + '/fram/' + nameextn
            pathanimextn = gdat.pathimag + rtag + '/' + strgpdfn + '/anim/' + nameextn
        
            try:
                listfile = fnmatch.filter(os.listdir(pathframextn), '*_swep*.pdf')
            except:
                print '%s failed.' % pathframextn
                continue
    
            listfiletemp = []
            for thisfile in listfile:
                listfiletemp.extend((thisfile.split('_')[0]).rsplit('/', 1))
            
            listname = list(set(listfiletemp))
            if len(listname) == 0:
                continue
            
            shuffle(listname)
    
            for name in listname:
                
                strgtemp = '%s*_swep*.pdf' % name
                listfile = fnmatch.filter(os.listdir(pathframextn), strgtemp)
                numbfile = len(listfile)
                liststrgextn = []
                for k in range(numbfile):
                    liststrgextn.append((listfile[k].split(name)[1]).split('_')[0])
                
                liststrgextn = list(set(liststrgextn))
                
                for k in range(len(liststrgextn)):
            
                    listfile = fnmatch.filter(os.listdir(pathframextn), name + liststrgextn[k] + '_swep*.pdf')
                    numbfile = len(listfile)
                    
                    indxfilelowr = 0
                    
                    if indxfilelowr < numbfile:
                        indxfileanim = np.arange(indxfilelowr, numbfile)
                    else:
                        continue
                        
                    indxfileanim = choice(indxfileanim, replace=False, size=indxfileanim.size)
                    
                    cmnd = 'convert -delay 20 -density 300 -quality 100 '
                    for n in range(indxfileanim.size):
                        cmnd += '%s%s ' % (pathframextn, listfile[indxfileanim[n]])
    
                    namegiff = '%s%s.gif' % (pathanimextn, name + liststrgextn[k])
                    cmnd += ' ' + namegiff
                    print 'Processing %s' % namegiff
                    if not os.path.exists(namegiff):
                        print 'Run: %s, pdf: %s' % (rtag, strgpdfn)
                        print 'Making %s animation...' % name
                        os.system(cmnd)
                    else:
                        print 'GIF already exists.'
                        pass
                    print
    
    pathoutprtag = retr_pathoutprtag(rtag)
    filestat = open(pathoutprtag + 'stat.txt', 'a')
    filestat.write('animfinl written.\n')
    filestat.close()
    

def plot_samp(gdat, gdatmodi, strgstat, strgmodl, strgphas, strgpdfn='post', gdatmock=None, booltile=False):
    
    print 'plot_samp()'

    backtype = getattr(gdat, strgmodl + 'backtype')
    boolbfun = getattr(gdat, strgmodl + 'boolbfun')
    lablelemextn = getattr(gdat, strgmodl + 'lablelemextn')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    
    strgpfix = retr_strgpfix(strgstat, strgmodl, strgpdfn=strgpdfn)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl)
    
    if not booltile:
        samp = getattr(gdatobjt, strgpfix + 'samp')
    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    numbdeflsubhplot = getattr(gdat, strgmodl + 'numbdeflsubhplot')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    numbdeflsingplot = getattr(gdat, strgmodl + 'numbdeflsingplot')
    if numbtrap > 0:
        boolelemlens = getattr(gdat, strgmodl + 'boolelemlens')
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
        liststrgfeatcorr = getattr(gdat, strgmodl + 'liststrgfeatcorr')
        spectype = getattr(gdat, strgmodl + 'spectype')
        elemtype = getattr(gdat, strgmodl + 'elemtype')
        calcerrr = getattr(gdat, strgmodl + 'calcerrr')
        boolelemsbrtextsbgrdanyy = getattr(gdat, strgmodl + 'boolelemsbrtextsbgrdanyy')
        boolelemdeflsubhanyy = getattr(gdat, strgmodl + 'boolelemdeflsubhanyy')
        boolelempsfnanyy = getattr(gdat, strgmodl + 'boolelempsfnanyy')
        liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
        if strgstat != 'pdfn':
            indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
            numbelem = [[] for l in indxpopl]
            for l in indxpopl:
                numbelem[l] = samp[indxfixpnumbelem[l]].astype(int)
    if strgstat != 'pdfn' and lensmodltype != 'none' and (strgmodl == 'fitt' and gdat.datatype == 'mock'):
        numbsingcomm = getattr(gdatobjt, strgpfix + 'numbsingcomm')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    if lensmodltype != 'none' or hostemistype != 'none':
        indxsersfgrd = getattr(gdat, strgmodl + 'indxsersfgrd')
    numbback = getattr(gdat, strgmodl + 'numbback')
    indxback = getattr(gdat, strgmodl + 'indxback')
    convdiff = getattr(gdat, strgmodl + 'convdiff')
    listnamediff = getattr(gdat, strgmodl + 'listnamediff')
    listnameecomtotl = getattr(gdat, strgmodl + 'listnameecomtotl')
    unifback = getattr(gdat, strgmodl + 'unifback')
    listnameback = getattr(gdat, strgmodl + 'listnameback')
    namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
    
    if gdatmodi is not None:
        strgswep = '_%09d' % gdatmodi.cntrswep
    else:
        strgswep = ''
    
    if not booltile:
        # data count maps
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpdata', i, m)
            ## residual count maps
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpresi', i, m)
        
        if gdat.numbpixl > 1:
            if numbtrap > 0:
                if boolelemlens:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelem', tdim=True)
        
        # temp -- restrict other plots to indxmodlelemcomp
        if gdat.enerbins:
            for specconvunit in gdat.listspecconvunit:
                if not boolbfun:
                    plot_sbrt(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, specconvunit)
    
    if numbtrap > 0:
        # element feature histograms
        if not (strgmodl == 'true' and gdat.datatype == 'inpt'):
            
            limtydat = gdat.limtydathistfeat

            for l in indxpopl:
                strgindxydat = 'pop%d' % l
                for strgfeat in liststrgfeatodim[l]:
                    if not (strgfeat == 'flux' or strgfeat == 'mcut' or strgfeat == 'deltllik' or strgfeat == 'defs' or strgfeat == 'nobj') and \
                                                                                (gdat.shrtfram and strgstat == 'this' and strgmodl == 'fitt'):
                        continue
                    indxydat = [l, slice(None)]
                    scalxdat = getattr(gdat, 'scal' + strgfeat + 'plot')
                    factxdat = getattr(gdat, 'fact' + strgfeat + 'plot')
                    lablxdat = getattr(gdat, 'labl' + strgfeat + 'totl')
                    limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxdat, getattr(gdat, 'maxm' + strgfeat) * factxdat]
                    
                    if gdat.numbpixl > 1:
                        listydattype = ['totl', 'sden']
                    else:
                        listydattype = ['totl']
                    for ydattype in listydattype:
                        
    
                        ## plot the surface density of elements
                        if ydattype == 'sden':
                            
                            # plot the surface density of elements only for the amplitude feature
                            if strgfeat != namefeatampl: 
                                continue
                            
                            if gdat.sdenunit == 'degr':
                                factydat = (pi / 180.)**2 / (2. * gdat.maxmgang)**2
                                lablydat = r'$\Sigma_{%s}$ [deg$^{-2}$]' % lablelemextn[l]
                            if gdat.sdenunit == 'ster':
                                factydat = 1. / (2. * gdat.maxmgang)**2
                                lablydat = r'$\Sigma_{%s}$ [sr$^{-2}$]' % lablelemextn[l]
                        
                        ## plot the total number of elements
                        if ydattype == 'totl':
                            factydat = 1.
                            lablydat = r'$N_{%s}$' % lablelemextn[l]
                    
                        if ydattype == 'totl' and not (gdat.datatype == 'inpt' and gdat.rtagmock is None):
                            listhisttype = ['hist', 'histptfn']
                        else:
                            listhisttype = ['hist']
                        
                        boolhistprio = not booltile
                        for histtype in listhisttype:
                            
                            if histtype == 'histptfn':
                                
                                if gdat.truenumbtrap == 0 or gdat.priofactdoff == 0.:
                                    continue

                                if strgfeat == 'specplot' or strgfeat == 'spec' or strgfeat == 'deflprof':
                                    continue
                            
                                if not strgfeat in gdat.trueliststrgfeat[l]:
                                    continue
                            
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'hist' + strgfeat + 'pop%d' % l, \
                                              'mean' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                              lablydat=lablydat, factxdat=factxdat, histodim=True, factydat=factydat, ydattype=ydattype, \
                                              scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, boolhistprio=boolhistprio, \
                                              #indxydat=indxydat, strgindxydat=strgindxydat, \
                                              nameinte='histodim/', histtype=histtype)
    
    if not booltile:
        if numbtrap > 0:
            # element feature correlations
            for l in gdat.fittindxpopl:
                if strgmodl != 'true' and gdat.refrinfo and gdat.boolasscrefr:
                    for strgfeat in gdat.fittliststrgfeatodim[l]:
                        if not (strgfeat == 'flux' or strgfeat == 'mass' or strgfeat == 'deltllik' or strgfeat == 'nobj') and \
                                                                                    (gdat.shrtfram and strgstat == 'this' and strgmodl == 'fitt'):
                            continue
                        for q in gdat.indxrefr:
                            if not l in gdat.refrindxpoplassc[q]:
                                continue
                            if gdat.refrnumbelem[q] == 0:
                                continue
                            if not strgfeat in gdat.refrliststrgfeat[q] or strgfeat in gdat.refrliststrgfeatonly[q][l]:
                                continue
                            plot_scatassc(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, q, l, strgfeat)
                            plot_scatassc(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, q, l, strgfeat, plotdiff=True)
                    
        if not (gdat.shrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            # plots
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    if numbpopl > 1:
                        if numbtrap > 0:
                            for l in indxpopl:
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpdata', i, m, indxpoplplot=l)
        
            ## histograms of the number of counts per pixel
            limtxdat = [gdat.minmcntpmodl, gdat.maxmcntpmodl]
            for m in gdat.indxevtt: 
                for nameecom in listnameecomtotl:
                    if gdat.numbpixl > 1:
                        for i in gdat.indxener:
                            name = 'histcntp' + nameecom + 'en%02devt%d' % (i, m)
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meancntpdata', \
                                                                     scalydat='logt', scalxdat='logt', lablxdat=gdat.lablcnts, histodim=True, \
                                                                     lablydat='$N_{pix}$', limtydat=[0.5, gdat.numbpixl], limtxdat=limtxdat)
                    else:
                        name = 'histcntp' + nameecom + 'evt%d' % m
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                                       name, 'meancntpdata', scalydat='logt', scalxdat='logt', lablxdat=gdat.lablcnts, histodim=True, \
                                                         lablydat='$N_{pix}$', limtydat=[0.5, gdat.numbener], limtxdat=limtxdat)

            
            for q in gdat.indxrefr:
                for l in gdat.fittindxpopl:
                    setattr(gdat, 'minmcmplpop%dpop%d' % (l, q), gdat.minmcmpl)
                    setattr(gdat, 'maxmcmplpop%dpop%d' % (l, q), gdat.maxmcmpl)
                    setattr(gdat, 'corrcmplpop%dpop%d' % (l, q), None)
                    setattr(gdat, 'factcmplpop%dpop%dplot' % (l, q), 1.)
                    setattr(gdat, 'scalcmplpop%dpop%d' % (l, q), 'self')
                    setattr(gdat, 'lablcmplpop%dpop%d' % (l, q), '$c_{%d%d}$' % (l, q))
                    
                    setattr(gdat, 'minmfdispop%dpop%d' % (q, l), gdat.minmfdis)
                    setattr(gdat, 'maxmfdispop%dpop%d' % (q, l), gdat.maxmfdis)
                    setattr(gdat, 'corrfdispop%dpop%d' % (q, l), None)
                    setattr(gdat, 'factfdispop%dpop%dplot' % (q, l), 1.)
                    setattr(gdat, 'scalfdispop%dpop%d' % (q, l), 'self')
                    setattr(gdat, 'lablfdispop%dpop%d' % (q, l), '$f_{%d%d}$' % (q, l))
                
            
            ## highest amplitude element
            # temp
            if numbtrap > 0:
                # completeness and false discovery rate
                if strgmodl != 'true' and gdat.boolasscrefr:
                    for strgclas in ['cmpl', 'fdis']:
                        nameinte = strgclas + 'odim/'
                        limtydat = [getattr(gdat, 'minm' + strgclas), getattr(gdat, 'maxm' + strgclas)]
                        for l in gdat.fittindxpopl:
                            for q in gdat.indxrefr:
                                if not l in gdat.refrindxpoplassc[q]:
                                    continue
                                if gdat.refrnumbelem[q] == 0 and strgclas == 'cmpl' or gdat.fittnumbtrap == 0 and strgclas == 'fdis':
                                    continue
                                if strgclas == 'cmpl':
                                    lablydat = getattr(gdat, 'labl' + strgclas + 'pop%dpop%d' % (l, q))
                                    strgindxydat = 'pop%dpop%d' % (l, q)
                                else:
                                    lablydat = getattr(gdat, 'labl' + strgclas + 'pop%dpop%d' % (q, l))
                                    strgindxydat = 'pop%dpop%d' % (q, l)
                                for strgfeat in gdat.refrliststrgfeat[q]:
                                    if strgfeat == 'etag':
                                        continue
                                    if strgclas == 'fdis' and not strgfeat in liststrgfeatodim[l]:
                                        continue
                                    if not strgfeat.startswith('spec') and not strgfeat.startswith('defl') \
                                                         and not strgfeat in gdat.refrliststrgfeatonly[q][l] and \
                                                         not (gdat.datatype == 'mock' and (strgfeat.endswith('pars') or strgfeat.endswith('nrel'))):
                                        factxdat = getattr(gdat, 'fact' + strgfeat + 'plot')
                                        lablxdat = getattr(gdat, 'labl' + strgfeat + 'totl')
                                        scalxdat = getattr(gdat, 'scal' + strgfeat + 'plot')
                                        limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxdat, getattr(gdat, 'maxm' + strgfeat) * factxdat]
                                        
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgclas + strgfeat + strgindxydat, \
                                                  'mean' + strgfeat, lablxdat=lablxdat, \
                                                  lablydat=lablydat, factxdat=factxdat, \
                                                  #plottype='errr', \
                                                  scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                  omittrue=True, nameinte=nameinte)

            if numbtrap > 0:
                alph = 0.1
                if strgmodl == 'true':
                    pathtemp = gdat.pathinit
                else:
                    if strgstat == 'this':
                        pathtemp = gdat.pathplotrtag + strgpdfn + '/fram/'
                    elif strgstat == 'mlik':
                        pathtemp = gdat.pathplotrtag + strgpdfn + '/finl/'
                    elif strgstat == 'pdfn':
                        pathtemp = gdat.pathplotrtag + strgpdfn + '/finl/'
                colr = retr_colr(gdat, strgstat, strgmodl, indxpopl=None)
        
                # transdimensional element features projected onto the data axes
                if not (strgstat == 'pdfn' and not gdat.boolcondcatl):
                    for l in indxpopl:
                        if elemtype[l] == 'lght':
                            # PS spectra
                            if strgstat == 'pdfn':
                                specplot = [np.empty((gdat.numbenerplot, gdat.numbstkscond))]
                                for r in gdat.indxstkscond:
                                    specplot[0][:, r] = gdat.dictglob['poststkscond'][r]['specplot'][0, :]
                            else:
                                specplot = getattr(gdatobjt, strgpfix + 'specplot')
                            
                            listxdat = []
                            listplottype = []
                            
                            for k in range(specplot[l].shape[-1]):
                                listxdat.append(gdat.meanenerplot)
                                listplottype.append('lghtline')
                            
                            for specconvunit in gdat.listspecconvunit:
                                listydat = []
                                
                                for k in range(specplot[l].shape[-1]):
                                    specplottemp = specplot[l]
                                    if strgmodl == 'true':
                                        specplottemp = np.copy(specplottemp[0, :, k])
                                    else:
                                        specplottemp = np.copy(specplottemp[:, k])
                                    if specconvunit[0] == 'en01':
                                        specplottemp *= gdat.meanenerplot
                                    if specconvunit[0] == 'en02':
                                        specplottemp *= gdat.meanenerplot**2
                                    if specconvunit[0] == 'en03':
                                        # temp
                                        pass
                                    listydat.append(specplottemp)
                                
                                lablydat = getattr(gdat, 'lablflux' + specconvunit[0] + specconvunit[1] + 'totl')
                                if specconvunit[1] == gdat.nameenerunit:
                                    factydat = 1.
                                else:
                                    factydat = getattr(gdat, 'fact' + specconvunit[1] + gdat.nameenerunit)
                                strgtemp = specconvunit[0] + specconvunit[1]
                                if specconvunit[0] == 'en03':
                                    strgtemp += specconvunit[2]
                                path = pathtemp + strgstat + 'specpop%d%s%s.pdf' % (l,  strgtemp, strgswep)
                                for ydat in listydat:
                                    ydat *= factydat
                                limtydat = [np.amin(gdat.minmspec) * factydat, np.amax(gdat.maxmspec) * factydat]
                                tdpy.util.plot_gene(path, listxdat, listydat, scalxdat='logt', scalydat='logt', \
                                                           lablxdat=gdat.lablenertotl, colr=colr, alph=alph, \
                                                           plottype=listplottype, limtxdat=[gdat.minmener, gdat.maxmener], lablydat=lablydat, \
                                                           limtydat=limtydat)
                    
                    if lensmodltype == 'elem' or lensmodltype == 'full':

                        ## deflection profiles
                        if gdat.variasca and gdat.variacut:
                            lablxdat = gdat.lablgangtotl
                            if strgstat == 'pdfn':
                                deflprof = [np.empty((gdat.numbanglfull, gdat.numbstkscond))]
                                asca = [np.empty(gdat.numbstkscond)]
                                acut = [np.empty(gdat.numbstkscond)]
                                for r in gdat.indxstkscond:
                                    deflprof[0][:, r] = gdat.dictglob['poststkscond'][r]['deflprof'][0, :]
                                    asca[0][r] = gdat.dictglob['poststkscond'][r]['asca'][0]
                                    acut[0][r] = gdat.dictglob['poststkscond'][r]['acut'][0]
                            else:
                                deflprof = getattr(gdatobjt, strgpfix + 'deflprof')
                                asca = getattr(gdatobjt, strgpfix + 'asca')
                                acut = getattr(gdatobjt, strgpfix + 'acut')

                            for l in range(len(deflprof)):
                                xdat = gdat.meananglfull * gdat.anglfact
                                listydat = []
                                listvlinfrst = []
                                listvlinseco = []
                                
                                if 'deflprof' in elemtype[l]:
                                    print 'l'
                                    print l
                                    print 'deflprof'
                                    print deflprof

                                    if strgmodl == 'true':
                                        deflproftemp = deflprof[l][0, :, :]
                                    else:
                                        deflproftemp = deflprof[l]
                                    
                                    for k in range(deflprof[l].shape[-1]):
                                        listydat.append(deflproftemp[:, k] * gdat.anglfact)
                                        if strgmodl == 'true':
                                            ascatemp = asca[l][0, k]
                                            acuttemp = acut[l][0, k]
                                        else:
                                            ascatemp = asca[l][k]
                                            acuttemp = acut[l][k]
                                        listvlinfrst.append(ascatemp * gdat.anglfact) 
                                        listvlinseco.append(acuttemp * gdat.anglfact)
                                    
                                    if lensmodltype == 'host':
                                        indxfixpbeinhost = getattr(gdat, strgmodl + 'indxfixpbeinhost')
                                        beinhost = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'samp', strgpdfn, indxvarb=indxfixpbeinhost)
                                        listydat.append(xdat * 0. + gdat.anglfact * beinhost)
                                    path = pathtemp + strgstat + 'deflsubhpop%d%s.pdf' % (l, strgswep)
                                    limtydat = [1e-3, 1.]
                                    limtxdat = [1e-3, 1.]
                                    tdpy.util.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', \
                                                                        lablxdat=lablxdat, drawdiag=True, limtydat=limtydat, \
                                                                        limtxdat=limtxdat, colr=colr, alph=alph, lablydat=r'$\alpha$ [$^{\prime\prime}$]', \
                                                                        listvlinfrst=listvlinfrst, listvlinseco=listvlinseco)
                        
                if gdat.datatype == 'mock':
                    # pulsar masses
                    lablxdat = gdat.lablgangtotl
                    factxdat = gdat.anglfact
                    for l in indxpopl:
                        if elemtype[l] == 'lghtpntspuls':
                            limtydat = [gdat.minmmassshel, gdat.maxmmassshel]
                            lablydat = gdat.lablmassshel
                            name = 'massshelpop%d' % l
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                                           lablxdat=lablxdat, lablydat=lablydat, factxdat=factxdat, limtydat=limtydat)

                    if lensmodltype != 'none':
                        ## radial mass budget
                        factxdat = gdat.anglfact
                        lablxdat = gdat.lablanglfromhosttotl
                        for namecalc in ['delt', 'intg']:
                            
                            # host mass
                            for e in indxsersfgrd:
                                strgisfr = 'isf%d' % e
                                limtydat = [gdat.minmmcut, getattr(gdat, 'maxmmasshost' + strgisfr + namecalc + 'bein')]
                                lablydat = getattr(gdat, 'lablmasshost' + strgisfr + namecalc + 'totl')
                                name = 'masshost%s%s' % (strgisfr, namecalc)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                                          lablxdat=lablxdat, lablydat=lablydat, factxdat=factxdat, limtydat=limtydat)
                            
                            if boolelemdeflsubhanyy:
                                # subhalo masses
                                limtydat = [gdat.minmmcut, getattr(gdat, 'maxmmasssubh' + namecalc + 'bein')]
                                lablydat = getattr(gdat, 'lablmasssubh' + namecalc + 'totl')
                                name = 'masssubh%s' % (namecalc)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                                          lablxdat=lablxdat, lablydat=lablydat, factxdat=factxdat, limtydat=limtydat)

                                # subhalo mass fraction
                                limtydat = [1e-3, 0.1]
                                lablydat = getattr(gdat, 'lablfracsubh' + namecalc + 'totl')
                                name = 'fracsubh%s' % (namecalc)
                                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'meananglhalf', scalydat='logt', \
                                                            lablxdat=lablxdat, lablydat=lablydat, factxdat=factxdat, limtydat=limtydat)

                alph = 0.1

                if boolelempsfnanyy:
                    ## PSF radial profile
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            indxydat = [i, slice(None), m]
                            strgindxydat = 'en%02devt%d' % (i, m)
                            lablxdat = gdat.lablgangtotl
                            limtydat= np.array([1e-3, 1e3]) * gdat.anglfact**2
                            plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'psfn', \
                                                            'binsangl', indxydat=indxydat, strgindxydat=strgindxydat, scalydat='logt', \
                                                            factxdat=gdat.anglfact, lablxdat=lablxdat, lablydat=r'$\mathcal{P}$', limtydat=limtydat)
        
                # internally and externally corrected element feature histograms
                if gdat.datatype == 'inpt' and strgstat == 'pdfn' and gdat.rtagmock is not None:
                    limtydat = gdat.limtydathistfeat
                    factydat = 1.
                    for l in indxpopl:
                        strgindxydat = 'pop%d' % l
                        for strgfeat in liststrgfeatodim[l]:
                            if strgfeat.startswith('aerr') or strgfeat == 'specplot' or strgfeat == 'spec' or strgfeat == 'deflprof':
                                continue
                            scalxdat = getattr(gdat, 'scal' + strgfeat + 'plot')
                            factxdat = getattr(gdat, 'fact' + strgfeat + 'plot')
                            lablxdat = getattr(gdat, 'labl' + strgfeat + 'totl')
                            limtxdat = [getattr(gdat, 'minm' + strgfeat) * factxdat, getattr(gdat, 'maxm' + strgfeat) * factxdat]
                            lablydat = r'$N_{%s}$' % lablelemextn[l]
                            for namecorr in ['incr', 'excr']:
                                nameinte = namecorr + 'odim/'
                                for qq in gdatmock.indxrefr:
                                    if namecorr == 'excr':
                                        if not strgfeat in gdat.fittliststrgfeatextr[l]:
                                            continue
                                        q = gdat.listnamerefr.index(strgfeat[-4:])
                                        if getattr(gdat, 'crex' + strgfeat + 'pop%dpop%dpop%d' % (q, qq, l)) is None:
                                            continue
                                        name = namecorr + strgfeat + 'pop%dpop%dpop%d' % (q, qq, l)
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'mean' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                          lablydat=lablydat, factxdat=factxdat, histodim=True, factydat=factydat, ydattype='totl', \
                                                          scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                          nameinte=nameinte)
               
                                    else:
                                        if strgfeat in gdat.fittliststrgfeatextr[l]:
                                            continue
                                        name = namecorr + strgfeat + 'pop%dpop%d' % (qq, l)
                                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, name, 'mean' + strgfeat, scalydat='logt', lablxdat=lablxdat, \
                                                      lablydat=lablydat, factxdat=factxdat, histodim=True, factydat=factydat, ydattype='totl', \
                                                      scalxdat=scalxdat, limtydat=limtydat, limtxdat=limtxdat, \
                                                      nameinte=nameinte)


    if not (gdat.shrtfram and strgstat == 'this' and strgmodl == 'fitt'):
        if numbtrap > 0:
            # element feature correlations
            liststrgelemtdimvarb = getattr(gdat, 'liststrgelemtdimvarb' + strgphas)
            for strgelemtdimtype in gdat.liststrgelemtdimtype:
                for strgelemtdimvarb in liststrgelemtdimvarb:
                    if strgelemtdimvarb.startswith('cmpl'):
                        continue
                    for l0 in indxpopl:
                        for strgfrst in liststrgfeat[l0]:
                            
                            if strgfrst == 'spec' or strgfrst == 'specplot' or strgfrst == 'deflprof':
                                continue

                            for strgseco in liststrgfeat[l0]:
                                
                                if strgseco == 'spec' or strgseco == 'specplot' or strgseco == 'deflprof':
                                    continue
                                
                                if not checstrgfeat(strgfrst, strgseco):
                                    continue
                                    
                                if strgelemtdimvarb.startswith('hist'):
                                    
                                    strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%d' % l0
                                    plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                        l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
                                else:
                                    if booltile:
                                        continue

                                    if strgfrst.startswith('aerr') or strgseco.startswith('aerr'):
                                        continue
                                    if strgelemtdimvarb.startswith('fdis'):
                                        for q in gdat.indxrefr:
                                            strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%dpop%d' % (q, l0)
                                            plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                            l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
                                    elif strgelemtdimvarb.startswith('excr') or strgelemtdimvarb.startswith('incr'):
                                        for qq in gdatmock.indxrefr:
                                            if strgelemtdimvarb.startswith('excr'):
                                                for q in gdat.indxrefr:
                                                    if getattr(gdat, 'crex' + strgfrst + strgseco + 'pop%dpop%dpop%d' % (q, qq, l0)) is None:
                                                        continue
                                                    strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%dpop%dpop%d' % (q, qq, l0)
                                                    plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                                l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
                                            else:
                                                if strgfrst[-4:] in gdat.listnamerefr and strgseco[-4:] in gdat.listnamerefr:
                                                    continue
                                                strgtotl = strgelemtdimvarb + strgfrst + strgseco + 'pop%dpop%d' % (qq, l0)
                                                plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, \
                                                                                            l0, strgfrst, strgseco, strgtotl, strgpdfn=strgpdfn)
        
            if not (gdat.datatype == 'mock' and (gdat.truenumbelemtotl == 0 or gdat.truemaxmnumbelemtotl == 0)):
                for q0 in gdat.indxrefr:
                    if booltile:
                        continue
                    for l0 in indxpopl:
                        for refrstrgfrst in gdat.refrliststrgfeat[q0]:
                            if refrstrgfrst == 'spec' or refrstrgfrst == 'specplot' or refrstrgfrst == 'deflprof' or refrstrgfrst == 'etag':
                                continue
                            if refrstrgfrst in gdat.refrliststrgfeatonly[q0][l0]:
                                continue
                            for refrstrgseco in gdat.refrliststrgfeat[q0]:
                                if refrstrgseco in gdat.refrliststrgfeatonly[q0][l0]:
                                    continue
                                if refrstrgseco == 'spec' or refrstrgseco == 'specplot' or refrstrgseco == 'deflprof' or refrstrgseco == 'etag':
                                    continue
                                
                                if not checstrgfeat(refrstrgfrst, refrstrgseco):
                                    continue
                                        
                                if refrstrgfrst.startswith('aerr') or refrstrgseco.startswith('aerr') or refrstrgfrst == 'specplot' or refrstrgseco == 'specplot':
                                    continue
                                
                                strgtotl = 'cmpl' + refrstrgfrst + refrstrgseco + 'pop%dpop%d' % (l0, q0)
                                
                                plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, 'bind', 'cmpl', \
                                                                                q0, refrstrgfrst, refrstrgseco, strgtotl, strgpdfn=strgpdfn)
            
    if not booltile:
        if not (gdat.shrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            # data and model count scatter
            for m in gdat.indxevttplot:
                if gdat.numbpixl > 1:
                    for i in gdat.indxener:
                        plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, m, indxenerplot=i)
                else:
                    plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, m)

            ## spatial priors
            # temp
            if gdat.numbpixl > 1:
                if numbtrap > 0:
                    liststrgfeatmodu = getattr(gdat, strgmodl + 'liststrgfeatmodu')
                    liststrgpdfnmodu = getattr(gdat, strgmodl + 'liststrgpdfnmodu')
                    for l in indxpopl:
                        for strgfeat, strgpdfn in zip(liststrgfeatmodu[l], liststrgpdfnmodu[l]):
                            if strgpdfn == 'tmplreln':
                                plot_genemaps(gdat, gdatmodi, 'fitt', strgpdfn, 'lpdfspatpriointp', tdim=True)
                            if strgpdfn == 'tmplgaum':
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'lpdfspatpriointp', tdim=True)
            
            # model count maps
            ## backgrounds
            if gdat.numbpixl > 1:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        for c in indxback:
                            if boolbfun:
                                continue
                            if not unifback[c]:
                                plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpback%04d' % c, i, m, strgcbar='cntpdata')
                
                ## count error
                if strgmodl != 'true':
                    if numbtrap > 0:
                        for l in indxpopl:
                            if calcerrr[l]:
                                for i in gdat.indxener:
                                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntperrr', i, -1, strgcbar='cntpresi')
                
                ## diffuse components 
                for i in gdat.indxener:
                    for k, name in enumerate(listnamediff):
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntp%s' % (name), i, strgcbar='cntpdata')
            
                ## model count maps
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpmodl', i, m, strgcbar='cntpdata')
            
                # likelihood
                if strgmodl != 'true':
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'llik', i, m, strgcbar='llikmaps')
                
                if lensmodltype != 'none':
                    ## lensing signal to noise
                    if strgmodl == 'true':
                        for i in gdat.indxener:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 's2nr', i, -1)
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magn', tdim=True)
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'conv', tdim=True)
                    for i in gdat.indxener:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplens', i, strgcbar='cntpdata', tdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplensgradmgtd', i, strgcbar='cntpdata', tdim=True)
            
            if gdat.penalpridiff:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                                'psecodimdatapntsen%02devt%d' % (i, m), 'meanmpolodim', lablxdat='$l$', lablydat='$P_{resi}(l)$', \
                                                                                                                 limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                        plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'psecodimdatapntsprioen%02devt%d' % (i, m), 'meanmpolodim', lablxdat='$l$', \
                                                                                           lablydat='$P_{prio}(l)$', limtydat=[1e-2, 2.], scalxdat='logt', scalydat='logt')
                
            if lensmodltype != 'none':
                indxydat = [slice(None)]
                strgindxydat = ''
                factxdat = 1e-3
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convpsecodim', 'meanwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P(k)$', limtydat=[1e-1, 1e2], \
                                                                          factxdat=factxdat, scalxdat='logt', scalydat='logt', indxydat=indxydat, strgindxydat=strgindxydat)
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'histdefl', 'meandefl', \
                                                                        scal='self', lablxdat=r'$\alpha$ [arcsec]', lablydat=r'$N_{pix}$', factxdat=gdat.anglfact, \
                                                                                 strgindxydat=strgindxydat, indxydat=indxydat, histodim=True)
            if numbtrap > 0 and boolelemdeflsubhanyy:
                indxydat = [slice(None)]
                strgindxydat = ''
                factxdat = 1e-3
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convpsecelemodim', 'meanwvecodim', lablxdat='$k$ [1/kpc]', lablydat='$P_{sub}(k)$', \
                                       factxdat=factxdat, strgindxydat=strgindxydat, indxydat=indxydat, limtydat=[1e-5, 1e-1], scalxdat='logt', scalydat='logt')
                plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'histdeflsubh', 'meandeflsubh', scal='self', lablxdat=r'$\alpha$ [arcsec]', \
                                       strgindxydat=strgindxydat, indxydat=indxydat, lablydat=r'$N_{pix}$', factxdat=gdat.anglfact, histodim=True)
            
            if lensmodltype != 'none':
                for i in gdat.indxener:
                    plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrd', i, -1, strgcbar='cntpdata')
                    if numbtrap > 0 and boolelemsbrtextsbgrdanyy:
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrdgalx', i, -1, strgcbar='cntpdata')
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntpbgrdexts', i, -1, strgcbar='cntpdata')
                
                # gradient of the lens emission
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'cntplensgrad', indxenerplot=i, indxevttplot=m)
                
        if True:
        # temp
        #if not (gdat.shrtfram and strgstat == 'this' and strgmodl == 'fitt'):
            if lensmodltype != 'none':
                # overall deflection field
                plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, multfact=0.1)
                
                # deflection field due to individual lenses
                for k in range(numbdeflsingplot):  
                    if k == 0:
                        multfact = 0.1
                    elif k == 1:
                        multfact = 1.
                    elif k >= 2:
                        multfact = 10.
                    plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, indxdefl=k, multfact=multfact)
                
                # residual deflection field
                if strgmodl == 'fitt' and gdat.datatype == 'mock':
                    plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgcomp='resi', multfact=100.)
                    if strgstat != 'pdfn':
                        for k in range(numbsingcomm):
                            plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgcomp='resi', indxdefl=k, multfact=100.)
                    
                    if gdat.numbpixl > 1:
                        if numbtrap > 0:
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelemresi', tdim=True)
                            plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'convelemresiperc', tdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magnresi', tdim=True)
                        plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, 'magnresiperc', tdim=True)
    

def dele_rtag(rtag):
    
    pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'
    pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
    
    cmnd = 'rm -rf %s%s' % (pathdata, rtag)
    print cmnd
    os.system(cmnd)
    cmnd = 'rm -rf %s%s' % (pathimag, rtag)
    os.system(cmnd)
    print cmnd


def plot_infopvks(gdat, gdatprio, name, namefull, nameseco=None):
    
    pvks = getattr(gdat, 'pvks' + namefull)

    info = getattr(gdat, 'info' + namefull)

    path = gdat.pathinfo + 'info' + namefull

    if nameseco is not None:
       
        limtfrst = getattr(gdat, 'limt' + name)
        limtseco = getattr(gdat, 'limt' + nameseco)
        factplotfrst = getattr(gdat, 'fact' + name + 'plot')
        factplotseco = getattr(gdat, 'fact' + nameseco + 'plot')
        varbfrst = getattr(gdat, 'mean' + name) * factplotfrst
        varbseco = getattr(gdat, 'mean' + nameseco) * factplotseco
        indxpoplfrst = int(namefull[-1])
        
        # information gain
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.pcolor(varbfrst, varbseco, info, cmap='Greys')
        plt.colorbar(imag)
        plot_sigmcont(gdat, 'fitt', axis, name, indxpoplfrst, strgseco=nameseco)
        scalfrst = getattr(gdat, 'scal' + name + 'plot')
        scalseco = getattr(gdat, 'scal' + nameseco + 'plot')
        if scalfrst == 'logt':
            axis.set_xscale('log')
        if scalseco == 'logt':
            axis.set_yscale('log')
        axis.set_xlabel(getattr(gdat, 'labl' + name + 'totl'))
        axis.set_ylabel(getattr(gdat, 'labl' + nameseco + 'totl'))
        axis.set_xlim(np.array(limtfrst) * factplotfrst)
        axis.set_ylim(np.array(limtseco) * factplotseco)
        plt.tight_layout()
        plt.savefig(path)
        plt.close(figr)

        # KS test p value
        pathpvkstdim = gdat.pathinfo + 'pvks' + namefull
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.pcolor(varbfrst, varbseco, pvks, cmap='Greys')
        plt.colorbar(imag)
        plot_sigmcont(gdat, 'fitt', axis, name, indxpoplfrst, strgseco=nameseco)
        if scalfrst == 'logt':
            axis.set_xscale('log')
        if scalseco == 'logt':
            axis.set_yscale('log')
        axis.set_xlabel(getattr(gdat, 'labl' + name + 'totl'))
        axis.set_ylabel(getattr(gdat, 'labl' + nameseco + 'totl'))
        axis.set_xlim(np.array(limtfrst) * factplotfrst)
        axis.set_ylim(np.array(limtseco) * factplotseco)
        plt.tight_layout()
        plt.savefig(pathpvkstdim)
        plt.close(figr)

    elif name != namefull:
        
        scal = getattr(gdat, 'scal' + name + 'plot')
        
        lablydat = '$D_{KL}$'
        lablxdat = getattr(gdat, 'labl' + name + 'totl')
        xdat = getattr(gdat, 'mean' + name) * getattr(gdat, 'fact' + name + 'plot')
        ydat = getattr(gdat, 'info' + namefull)
        tdpy.mcmc.plot_plot(path, xdat, ydat, lablxdat, lablydat, scal)
        
        ydat = getattr(gdat, 'pvks' + namefull)
        pathpvks = gdat.pathinfo + 'pvks' + namefull
        tdpy.mcmc.plot_plot(pathpvks, xdat, ydat, lablxdat, '$p_{KS}$', scal)
        
    else:
        # horizontal axis
        xdat = getattr(gdat, 'mean' + name) * getattr(gdat, 'fact' + name + 'plot')
        lablxdat = getattr(gdat, 'labl' + name + 'totl')
        
        # scaling
        scal = getattr(gdat, 'scal' + name) 
        
        # common title
        titl = '$D_{KL} = %.3g$, KS = %.3g $\sigma$' % (info, pvks)

        # DKL density
        pathdinf = gdat.pathinfo + 'dinf' + namefull
        ydat = getattr(gdat, 'infodens' + namefull)
        lablydat = r'$\rho_{D_{KL}}$'
        tdpy.mcmc.plot_plot(pathdinf, xdat, ydat, lablxdat, lablydat, scal, titl=titl)
        
        # prior and posterior PDFs
        pathpdfn = gdat.pathinfo + 'pdfn' + namefull
        lablydat = r'$P$'
        ydat = [getattr(gdat, 'pdfnpost' + namefull), getattr(gdatprio, 'pdfnprio' + namefull)]
        legd = ['$P$(%s|$D$)' % lablxdat, '$P$(%s)' % lablxdat]
        tdpy.mcmc.plot_plot(pathpdfn, xdat, ydat, lablxdat, lablydat, scal, colr=['k', 'k'], linestyl=['-', '--'], legd=legd, titl=titl)


def plot_finl(gdat=None, gdatprio=None, rtag=None, strgpdfn='post', gdatmock=None, booltile=None):
    
    if gdat.verbtype > 0:
        print 'plot_finl()'
        print 'Producing postprocessing plots...'

    timetotlinit = gdat.functime()
    
    gdat.strgbest = 'ML'
    
    if not booltile:
        # terms in the log-acceptance probability
        listindxsamptotlproptotl = getattr(gdat, 'list' + strgpdfn + 'indxsamptotlproptotl')
        listindxsamptotlpropaccp = getattr(gdat, 'list' + strgpdfn + 'indxsamptotlpropaccp')
        listindxsamptotlpropreje = getattr(gdat, 'list' + strgpdfn + 'indxsamptotlpropreje')
        for n in gdat.indxproptype:
            pathbase = getattr(gdat, 'path' + strgpdfn + 'finl%s' % gdat.nameproptype[n])
            for k in gdat.indxtermlacp:
                varb = getattr(gdat, 'list' + strgpdfn + gdat.listnametermlacp[k])
                labl = gdat.listlabltermlacp[k]
                
                if listindxsamptotlproptotl[n].size > 0 and (varb[listindxsamptotlproptotl[n]] != 0.).any():
                    path = pathbase + gdat.listnametermlacp[k] + 'totl'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlproptotl[n]], labl, titl=gdat.nameproptype[n] + ', Total')
                
                if listindxsamptotlpropaccp[n].size > 0 and (varb[listindxsamptotlpropaccp[n]] != 0.).any():
                    path = pathbase + gdat.listnametermlacp[k] + 'accp'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlpropaccp[n]], labl, titl=gdat.nameproptype[n] + ', Accepted')
                
                if listindxsamptotlpropreje[n].size > 0 and (varb[listindxsamptotlpropreje[n]] != 0.).any():
                    path = pathbase + gdat.listnametermlacp[k] + 'reje'
                    tdpy.mcmc.plot_trac(path, varb[listindxsamptotlpropreje[n]], labl, titl=gdat.nameproptype[n] + ', Rejected')
            
        if gdat.checprio and strgpdfn == 'post' and not booltile:
            # this works only for scalar variables -- needs to be generalized to all variables
            if gdatprio is None:
                pathoutprtag = retr_pathoutprtag(rtag)
                path = pathoutprtag + 'gdatfinlprio'
                gdatprio = readfile(path)

            for namevarbscal in gdat.listnamevarbscal:
                plot_infopvks(gdat, gdatprio, namevarbscal, namevarbscal)
            for l in gdat.fittindxpopl:
                for strgfeatfrst in gdat.fittliststrgfeat[l]:
                    if strgfeatfrst == 'spec' or strgfeatfrst == 'deflprof' or strgfeatfrst == 'specplot':
                        continue
                    plot_infopvks(gdat, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + 'pop%d' % l)
                    for strgfeatseco in gdat.fittliststrgfeat[l]:
                        if strgfeatseco == 'spec' or strgfeatseco == 'deflprof' or strgfeatseco == 'specplot':
                            continue
                        
                        if not checstrgfeat(strgfeatfrst, strgfeatseco):
                            continue
                                        
                        plot_infopvks(gdat, gdatprio, strgfeatfrst, 'hist' + strgfeatfrst + strgfeatseco + 'pop%d' % l, nameseco=strgfeatseco)
        
        listsamp = getattr(gdat, 'list' + strgpdfn + 'samp')
        listsamp = getattr(gdat, 'list' + strgpdfn + 'samp')
        listfixp = getattr(gdat, 'list' + strgpdfn + 'fixp')
    
        listboolpropfilt = getattr(gdat, 'list' + strgpdfn + 'boolpropfilt')
        listmemoresi = getattr(gdat, 'list' + strgpdfn + 'memoresi')
        listindxproptype = getattr(gdat, 'list' + strgpdfn + 'indxproptype')
        listsampproc = getattr(gdat, 'list' + strgpdfn + 'sampproc')
    
        # Gelman-Rubin test
        pathdiag = getattr(gdat, 'path' + strgpdfn + 'finldiag')
        if gdat.numbproc > 1:
            if np.isfinite(gdat.gmrbstat).all():
                if gdat.verbtype > 0:
                    print 'Gelman-Rubin TS...'
        
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                minm = min(np.amin(gdat.gmrbstat), np.amin(gdat.gmrbfixp))
                maxm = max(np.amax(gdat.gmrbstat), np.amax(gdat.gmrbfixp))
                bins = np.linspace(minm, maxm, 40)
                axis.hist(gdat.gmrbstat.flatten(), bins=bins, label='Data proj.')
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
                        maps = gdat.gmrbstat[i, :, m]
                        path = pathdiag + 'gmrbdataen%02devt%d.pdf' % (i, m)
                        tdpy.util.plot_maps(path, maps, indxpixlrofi=gdat.indxpixlrofi, numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                   minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                   minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
            else:
                print 'Inappropriate Gelman-Rubin test statistics encountered.'
    
        # plot autocorrelation
        if gdat.verbtype > 0:
            print 'Autocorrelation...'
        tdpy.mcmc.plot_atcr(pathdiag, gdat.atcrcntp[0, 0, 0, 0, :], gdat.timeatcrcntp[0, 0, 0, 0], strgextn='cntp')
        tdpy.mcmc.plot_atcr(pathdiag, gdat.atcrpara[0, 0, :], gdat.timeatcrpara[0, 0], strgextn='para')
        print 'Autocorrelation times:'
        for k, namepara in enumerate(gdat.fittnamepara):
            print '%s %g' % (namepara, np.mean(gdat.timeatcrpara[:, k]))
        
        # plot proposal efficiency
        if gdat.verbtype > 0:
            print 'Acceptance ratio...'
        numbtimemcmc = 20
        binstimemcmc = np.linspace(0., gdat.numbswep, numbtimemcmc)
        numbtick = 2
        sizefigrydat = 0.75 * gdat.numbproptype * gdat.plotsize
        figr, axgr = plt.subplots(gdat.numbproptype, 1, figsize=(gdat.plotsize, sizefigrydat), sharex='all')
        if gdat.numbproptype == 1:
            axgr = [axgr]
        for k, axis in enumerate(axgr):
            histtotl = axis.hist(listindxsamptotlproptotl[n+k], bins=binstimemcmc)[0]
            histaccp = axis.hist(listindxsamptotlpropaccp[n+k], bins=binstimemcmc)[0]
            axis.set_ylabel('%s' % gdat.nameproptype[k])
            if k == gdat.numbproptype - 1:
                axis.set_xlabel('$i_{samp}$')
        plt.tight_layout()
        figr.savefig(pathdiag + 'accpratiproptype.pdf')
        plt.close(figr)
   
        if gdat.verbtype > 0:
            print 'Proposal execution times...'
        
        ## time performance
        #listchro = np.empty((gdat.numbswep, gdat.numbchro))
        #listchro = []
        #for k, name in enumerate(gdat.listnamechro):
        #    #listchro[:, k] = getattr(gdat, 'list' + strgpdfn + 'chro' + name).flatten() * 1e3
        #    listchro.append(getattr(gdat, 'list' + strgpdfn + 'chro' + name).flatten() * 1e3)
        #pathdiag = getattr(gdat, 'path' + strgpdfn + 'finldiag')
        #figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
        #axis.violin(listchro)
        #axis.set_yscale('log')
        #axis.set_ylabel('$t$ [ms]')
        #axis.set_xticklabels(gdat.listlablchro)
        #axis.axvline(mean(chro), ls='--', alpha=0.2, color='black')
        #figr.savefig(pathdiag + 'chro.pdf' % gdat.listnamechro[k])
        #plt.close(figr)

    # temp
    gdat.legdpmea = 'Mean'

    # posterior versions of the frame plots
    plot_samp(gdat, None, 'pdfn', 'fitt', 'finl', strgpdfn=strgpdfn, gdatmock=gdatmock, booltile=booltile)
   
    if booltile:
        return

    if gdat.fittnumbtrap > 0:
        if gdat.verbtype > 0:
            print 'A mosaic of samples...'
    
        ## mosaic of images of posterior catalogs
        if gdat.numbpixl > 1:
            plot_mosa(gdat, strgpdfn)
    
    ## randomly selected trandimensional parameters
    if gdat.fittnumbtrap > 0:
        if gdat.verbtype > 0:
            print 'Transdimensional parameters...'
    
        # choose the parameters based on persistence
        stdvlistsamptran = np.std(listsamp[:, gdat.fittindxsamptrap], axis=0)
        indxtrapgood = np.where(stdvlistsamptran > 0.)[0]
        numbtrapgood = indxtrapgood.size
        numbtrapplot = min(3, numbtrapgood)
        if numbtrapplot > 0:
            indxtrapplot = np.sort(choice(gdat.fittindxsamptrap[indxtrapgood], size=numbtrapplot, replace=False))

            path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
            tdpy.mcmc.plot_grid(path, 'listelemfrst', listsamp[:, gdat.fittindxsamptrap[:3]] * gdat.fittfactplotpara[gdat.fittindxsamptrap[:3]], \
                                                                                                [gdat.fittlablpara[k] for k in gdat.fittindxsamptrap[:3]])

            path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
            tdpy.mcmc.plot_grid(path, 'listsamp', listsamp[:, indxtrapplot], ['%d' % k for k in indxtrapplot])
            path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
            tdpy.mcmc.plot_grid(path, 'listsamp', listsamp[:, indxtrapplot] * gdat.fittfactplotpara[indxtrapplot], \
                                                                                            [gdat.fittlablpara[k] for k in indxtrapplot])
    
    if gdat.verbtype > 0:
        print 'Scalar variables...'
    # scalar variables
    ## trace and marginal distribution of each parameter
    for name in gdat.listnamevarbscal:
        
        if gdat.verbtype > 0:
            print 'Working on %s...' % name
        scal = getattr(gdat, 'scal' + name) 
        factplot = getattr(gdat, 'fact' + name + 'plot')
        corr = getattr(gdat, 'corr' + name)
        if corr is None:
            truepara = None
        else:
            truepara = getattr(gdat, 'corr' + name) * factplot
        labltotl = getattr(gdat, 'labl' + name + 'totl')
        
        listvarb = getattr(gdat, 'list' + strgpdfn + name) * factplot
        if listvarb.ndim != 1:
            if listvarb.shape[1] == 1:
                listvarb = listvarb[:, 0]
            else:
                raise Exception('')
        
        mlik = getattr(gdat, 'mlik' + name) * factplot
        path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscaltrac') + name
        tdpy.mcmc.plot_trac(path, listvarb, labltotl, truepara=truepara, scalpara=scal, varbdraw=[mlik], labldraw=[''], colrdraw=['r'])
        path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalhist') + name
        tdpy.mcmc.plot_hist(path, listvarb, labltotl, truepara=truepara, scalpara=scal, varbdraw=[mlik], labldraw=[''], colrdraw=['r'])
       
        for nameseco in gdat.listnamevarbscal:
            
            if name == nameseco:
                continue
            
            if gdat.verbtype > 0:
                print 'Working on correlation of %s with %s...' % (name, nameseco)
            
            pathjoin = getattr(gdat, 'path' + strgpdfn + 'finlvarbscaljoin')
            scalseco = getattr(gdat, 'scal' + nameseco) 
            factplotseco = getattr(gdat, 'fact' + nameseco + 'plot')
            corrseco = getattr(gdat, 'corr' + nameseco)
            if corrseco is None:
                trueparaseco = None
            else:
                trueparaseco = getattr(gdat, 'corr' + nameseco) * factplotseco
            labltotlseco = getattr(gdat, 'labl' + nameseco + 'totl')
            listvarbseco = getattr(gdat, 'list' + strgpdfn + nameseco) * factplotseco
            mlikseco = getattr(gdat, 'mlik' + nameseco) * factplotseco
            
            if listvarbseco.ndim != 1:
                if listvarbseco.shape[1] == 1:
                    listvarbseco = listvarbseco[:, 0]
                else:
                    raise Exception('')
                
            listjoin = np.vstack((listvarb, listvarbseco)).T
    
            print 'truepara'
            print truepara
            print 'scal'
            print scal
            print 'scalseco'
            print scalseco
            print 'mlik'
            print mlik
            print 'trueparaseco'
            print trueparaseco
            tdpy.mcmc.plot_grid(pathjoin, name + nameseco, listjoin, [labltotl, labltotlseco], scalpara=[scal, scalseco], truepara=[truepara, trueparaseco], \
                                                                                                                   join=True, varbdraw=[mlik, mlikseco])

    if gdat.verbtype > 0:
        print 'Fixed dimensional parameter covariance...'
    
    ### covariance
    ## overall
    path = getattr(gdat, 'path' + strgpdfn + 'finlvarbscalcova')
    truepara = gdat.fittcorrfixp * gdat.fittfactfixpplot
    mlikpara = gdat.mlikfixp * gdat.fittfactfixpplot
    tdpy.mcmc.plot_grid(path, 'fixp', listfixp * gdat.fittfactfixpplot[None, :], gdat.fittlablfixptotl, truepara=truepara, varbdraw=mlikpara)
    
    # stacked posteiors binned in position and flux
    if gdat.fittnumbtrap > 0 and gdat.numbpixl > 1:
        liststrgbins = ['quad', 'full']
        for l in gdat.fittindxpopl:
            plot_histlgalbgalelemstkd(gdat, strgpdfn, l, 'cumu')
            for strgbins in liststrgbins:
                for strgfeatsign in gdat.fittliststrgfeatsign[l]:
                    plot_histlgalbgalelemstkd(gdat, strgpdfn, l, strgbins, strgfeatsign)

    if gdat.verbtype > 0:
        print 'Prior and likelihood...'
    
    for strgpdfntemp in ['lpritotl', 'lliktotl']:

        if strgpdfntemp == 'lpritotl':
            labltemp = '\ln P(M)'
        if strgpdfntemp == 'lliktotl':
            labltemp = '\ln P(D|M)'
        labl = r'$%s$' % labltemp

        path = getattr(gdat, 'path' + strgpdfn + 'finl') + strgpdfntemp
        
        varb = getattr(gdat, 'list' + strgpdfn + strgpdfntemp)
        tdpy.mcmc.plot_hist(path, varb, labl)
        varbdraw = []
        labldraw = []
        colrdraw = []
        if gdat.datatype == 'mock':
            varbdraw += [getattr(gdat, 'true' + strgpdfntemp)]
            labldraw += ['True model']
            colrdraw += [gdat.refrcolr]
        
        tdpy.mcmc.plot_trac(path, getattr(gdat, 'list' + strgpdfn + strgpdfntemp), labl, varbdraw=varbdraw, labldraw=labldraw, colrdraw=colrdraw)
    
    # plot resident memory
    figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))
    axis.plot(gdat.indxswep, np.mean(listmemoresi, 1) / float(2**30))
    axis.set_ylabel(r'$M$ [GB]')
    axis.set_xlabel(r'$i_{samp}$')
    plt.tight_layout()
    figr.savefig(pathdiag + 'memoresi.pdf')
    plt.close(figr)

    timetotlfinl = gdat.functime()
    if gdat.verbtype > 0:
        print 'Plots and animations are produced in %.3g seconds.' % (timetotlfinl - timetotlinit)


def plot_sbrt(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, specconvunit):
    
    strgpfix = retr_strgpfix(strgstat, strgmodl, strgpdfn=strgpdfn)
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl)
    
    backtype = getattr(gdat, strgmodl + 'backtype')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    specback = getattr(gdat, strgmodl + 'specback')
    indxback = getattr(gdat, strgmodl + 'indxback')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    boolelemsbrtdfncanyy = getattr(gdat, strgmodl + 'boolelemsbrtdfncanyy')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    indxbacpback = getattr(gdat, strgmodl + 'indxbacpback')
    numblablsbrt = getattr(gdat, strgmodl + 'numblablsbrt')
    numblablsbrtspec = getattr(gdat, strgmodl + 'numblablsbrtspec')
    if numbtrap > 0:
        boolelemlghtanyy = getattr(gdat, strgmodl + 'boolelemlghtanyy')
    
    if lensmodltype != 'none' or hostemistype != 'none':
        indxsersfgrd = getattr(gdat, strgmodl + 'indxsersfgrd')
    samp = getattr(gdatobjt, strgpfix + 'samp')
    if numbtrap > 0:
        indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
    
    if gdat.numbpixl == 1:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    
    if specconvunit[1] == gdat.nameenerunit:
        factydat = 1.
    else:
        factydat = getattr(gdat, 'fact' + specconvunit[1] + gdat.nameenerunit)

    for b, namespatmean in enumerate(gdat.listnamespatmean):
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        
        # plot reference spectra
        if gdat.listprefsbrtlabl is not None:
            for k in range(len(gdat.listprefsbrtlabl)):
                if gdat.listprefsbrttype[k] == 'shad':
                    factenerrefr = [[] for a in range(3)]
                    for a in range(3):
                        factenerrefr[a] = retr_factener(specconvunit[0], gdat.listprefsbrtener[k][a])
                    axis.plot(gdat.listprefsbrtener[k][0], factydat * gdat.listprefsbrtsbrt[k][0] * factenerrefr[0], color='m', label=gdat.listprefsbrtlabl[k])
                    enerpoly = np.empty(gdat.listprefsbrtener[k][1].size + gdat.listprefsbrtener[k][2].size)
                    enerpoly[:gdat.listprefsbrtener[k][1].size] = gdat.listprefsbrtener[k][1]
                    enerpoly[gdat.listprefsbrtener[k][1].size:] = gdat.listprefsbrtener[k][2][::-1]
                    sbrtpoly = np.empty(gdat.listprefsbrtener[k][1].size + gdat.listprefsbrtener[k][2].size)
                    sbrtpoly[:gdat.listprefsbrtener[k][1].size] = factydat * gdat.listprefsbrtsbrt[k][1] * factenerrefr[1]
                    sbrtpoly[gdat.listprefsbrtener[k][1].size:] = factydat * gdat.listprefsbrtsbrt[k][2][::-1] * factenerrefr[2][::-1]
                    axis.fill(enerpoly, sbrtpoly, color='m', alpha=0.5)
                else:
                    factenerrefr = retr_factener(specconvunit[0], gdat.listprefsbrtener[k][1])
                    axis.errorbar(gdat.listprefsbrtener[k][1], factydat * gdat.listprefsbrtsbrt[k][1] * factenerrefr, label=gdat.listprefsbrtlabl[k], color='m')
        
        if strgmodl == 'true':
            liststrgmodl = [strgmodl]
            listgdatobjt = [gdat]
        if strgmodl == 'fitt' and (strgstat == 'this' or strgstat == 'pdfn'):
            if gdat.datatype == 'mock':
                liststrgmodl = [strgmodl, 'true']
                listgdatobjt = [gdatobjt, gdat]
            else:
                liststrgmodl = [strgmodl]
                listgdatobjt = [gdatobjt]
        numbstrgstattemp = len(liststrgmodl)
        for a in range(numbstrgstattemp):
            
            indxploteleminit = []
            indxplotelemendd = []
                
            listlegdsbrtspec = getattr(gdat, strgmodl + 'listlegdsbrtspec')

            # number of transdimensional elements to be overplotted
            numbelemtemp = 0
            
            if gdat.numbpixl == 1 and strgstat != 'pdfn':
                if liststrgmodl[a] == 'fitt':
                    numbelem = [[] for l in indxpopl]
                    for l in indxpopl:
                        numbelem[l] = samp[indxfixpnumbelem[l]].astype(int)
                        numbelemtemp += np.sum(numbelem[l])
                else:
                    for q in gdat.indxrefr:
                        numbelemtemp += np.sum(gdat.refrnumbelem[q])
                
            numbplot = numblablsbrtspec + numbelemtemp
            listydat = np.zeros((numbplot, gdat.numbener))
            listyerr = np.zeros((2, numbplot, gdat.numbener))
            
            cntr = 0
            cntrdata = cntr

            ## data
            listydat[cntr, :] = gdat.sbrtdatamean[b]
            listyerr[:, cntr, :] = gdat.sbrtdatastdv[b]
            cntr += 1
            
            for c in indxback:
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtback%04dmea%d' % (c, b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtback%04dmea%d' % (c, b), strgpdfn, strgmome='errr')
                cntr += 1
            
            if numbtrap > 0 and boolelemsbrtdfncanyy and not (liststrgmodl[a] == 'true' and gdat.refrnumbelemtotl == 0):
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
            
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncsubtmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtdfncsubtmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
            
            if hostemistype != 'none':
                for e in indxsersfgrd:
                    listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrthostisf%dmea%d' % (e, b), strgpdfn)
                    if strgstat == 'pdfn':
                        listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrthostisf%dmea%d' % (e, b), strgpdfn, strgmome='errr')
                    cntr += 1
            
            if lensmodltype != 'none':
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtlensmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtlensmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
            
            if gdat.numbpixl == 1 and strgstat != 'pdfn':
                cntrline = cntr
                indxploteleminit.append(cntr)
                for l in indxpopl:
                    if liststrgmodl[a] == 'true':
                        for k in range(gdat.truenumbelem[l]):
                            listydat[cntr, :] = getattr(listgdatobjt[a], liststrgmodl[a] + 'spec')[l][0, :, k]
                            
                            if cntr == cntrline:
                                listlegdsbrtspec = listlegdsbrtspec[:cntr] + ['Lines'] + listlegdsbrtspec[cntr:]
                            else:
                                listlegdsbrtspec = listlegdsbrtspec[:cntr] + [None] + listlegdsbrtspec[cntr:]
                            
                            cntr += 1
                            if k == gdat.truenumbelem[l] - 1:
                                indxplotelemendd.append(k)
                    else:   
                        for k in range(numbelem[l]):
                            listydat[cntr, :] = getattr(listgdatobjt[a], strgstat + 'spec')[l][:, k]
                            
                            if cntr == cntrline:
                                listlegdsbrtspec = listlegdsbrtspec[:cntr] + ['Lines'] + listlegdsbrtspec[cntr:]
                            else:
                                listlegdsbrtspec = listlegdsbrtspec[:cntr] + [None] + listlegdsbrtspec[cntr:]
                
                            cntr += 1
                            if k == numbelem[l] - 1:
                                indxplotelemendd.append(k)
            ## total model
            if numblablsbrt > 1:
                listydat[cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtmodlmea%d' % (b), strgpdfn)
                if strgstat == 'pdfn':
                    listyerr[:, cntr, :] = retr_fromgdat(gdat, gdatmodi, strgstat, liststrgmodl[a], 'sbrtmodlmea%d' % (b), strgpdfn, strgmome='errr')
                cntr += 1
           
            if liststrgmodl[a] == 'true':
                listyerr = np.zeros((2, numbplot, gdat.numbener))
            
            # plot energy spectra of the data, background model components and total background
            if gdat.numbener > 1:
                
                listmrkr = ['o', '>', 's', 'h', '*', 'p', 'x']
                for k in range(100):
                    listmrkr.append('x')

                # determine the energy scaling factor
                if specconvunit[0] == 'en00':
                    factener = 1.
                if specconvunit[0] == 'en01':
                    factener = gdat.meanener
                if specconvunit[0] == 'en02':
                    factener = gdat.meanener**2
                if specconvunit[0] == 'en03':
                    # temp
                    pass
                    factener = 1.
                    #indxenerintv = np.where((gdat.meanener < specconvunit[4]) & (gdat.meanener > specconvunit[3]))[0]
                    #ener = np.concatenate((np.array([specconvunit[3]]), gdat.meanener[indxenerintv], np.array([specconvunit[4]])))
                    #
                    #for k in range(3):
                    #    if k == 0:
                    #        ydattemp = 
                    #    ydatminmener = np.interp(specconvunit[3], gdat.meanener, ydat)
                    #    ydatmaxmener = np.interp(specconvunit[4], gdat.meanener, ydat)
                    #    ydat = np.concatenate((np.array([ydatminmener]), ydat[indxenerintv], np.array([ydatmaxmener])))
                    #    ydat = np.trapz(ydat, gdat.meanener)
                    #
                    #yerrminmener = np.interp(specconvunit[3], gdat.meanener, yerr, axis=1)
                    #yerrmaxmener = np.interp(specconvunit[4], gdat.meanener, yerr, axis=1)
                    #ydat = np.stack((np.array([yerrminmener]), ydat[indxenerintv], np.array([yerrmaxmener])))
                    #
                    #
                    #yerr = np.trapz(yerr, gdat.meanener)


                xdat = gdat.meanener
                cntr = 0
                
                for k in range(listydat.shape[0]):
                    mrkr = listmrkr[cntr]
                    if k == cntrdata:
                        colr = 'black'
                        alph = 1.
                        linestyl = '-'
                    else:
                        colr = retr_colr(gdat, strgstat, liststrgmodl[a], indxpopl=None)
                        linestyl = '--'
                        alph = 0.5
                   
                    ydat = np.copy(listydat[k, :])
                    yerr = np.copy(listyerr[:, k, :])
                    
                    ydat *= factener
                    yerr *= factener
                    
                    ydat *= factydat
                    yerr *= factydat
                    
                    if k == cntrdata and a > 0:
                        continue
                    
                    if liststrgmodl[a] == 'fitt':
                        legd = listlegdsbrtspec[k]
                    else:
                        legd = None
                    
                    temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, color=colr, marker=mrkr, ls=linestyl, markersize=10, alpha=alph, label=legd)
                    for caps in listcaps:
                        caps.set_markeredgewidth(1)

                    if gdat.numbpixl == 1 and strgstat != 'pdfn':
                        if cntr != cntrline or k in indxplotelemendd:
                            cntr += 1
                    else:
                        cntr += 1

        if gdat.numbener > 1:
            axis.set_xlim([np.amin(gdat.binsener), np.amax(gdat.binsener)])
            
            if gdat.exprtype == 'chan':
                factminm = 1e-1
                factmaxm = 1e2
            elif gdat.exprtype == 'ferm':
                factminm = 1e1
                factmaxm = 1e-1
            else:
                factminm = 1e-4
                factmaxm = 1e0
            minmydat = factminm * gdat.factylimsbrt[0] * np.amax(listydat[cntrdata, :] * factener) * factydat
            maxmydat = factmaxm * gdat.factylimsbrt[1] * np.amax(listydat[cntrdata, :] * factener) * factydat
            limtydat = [minmydat, maxmydat]
            axis.set_ylim(limtydat)
            axis.set_yscale('log')
            axis.set_xlabel(gdat.lablenertotl)
            axis.set_xscale('log')
            labl = getattr(gdat, 'lablsbrt' + specconvunit[0] + specconvunit[1] + 'stertotl')
            axis.set_ylabel(labl)
            make_legd(axis, numbcols=2)
            
            plt.tight_layout()
            path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, 'sdenmean%s%s%s' % (namespatmean, specconvunit[0], specconvunit[1]))
            figr.savefig(path)
            plt.close(figr)
        

def retr_factener(strgconvunit, ener):
    
    if strgconvunit == 'en00':
        factener = np.ones_like(ener)
    
    if strgconvunit == 'en01':
        factener = ener
    
    if strgconvunit == 'en02':
        factener = ener**2
    
    if strgconvunit == 'en03':
        # temp
        pass
        factener = np.ones_like(ener)
    
    return factener


def plot_pdfntotlflux():

    minm = 1e-9
    maxm = 10e-9
    numbvarb = 90
    numbsamp = 100000
    numbbins = 40
    alph = 0.5
    
    binssing = np.linspace(minm, maxm, numbvarb + 1)
    meansing = (binssing[:-1] + binssing[1:]) / 2.
    deltsing = binssing[1:] - binssing[:-1]
    
    binsdoub = np.linspace(2. * minm, 2. * maxm, 2 * numbvarb)
    meandoub = (binsdoub[:-1] + binsdoub[1:]) / 2.
    deltdoub = binsdoub[1:] - binsdoub[:-1]
    
    bins = np.linspace(minm, 2. * maxm, 2 * numbvarb + 1)
    
    arry = np.empty((2, numbsamp))
    
    minmslop = 1.5
    maxmslop = 3.
    numbslop = 4
    sloparry = np.linspace(minmslop, maxmslop, numbslop)
    for n in range(numbslop):
        slop = sloparry[n]
        for k in range(2):
            arry[k, :] = (rand(numbsamp) * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))
        
        totl = np.sum(arry, 0)
        
        powrprob = (1. - slop) / (maxm**(1. - slop) - minm**(1. - slop)) * meansing**(-slop)
        
        convprob = convolve(powrprob, powrprob) * deltdoub[0]
        
        indxdoub = np.where(meandoub <= maxm)[0]
        convprobpoly = polyval(polyfit(meandoub[indxdoub], convprob[indxdoub], 8), meandoub[indxdoub])
        
        figr, axis = plt.subplots()
        axis.hist(arry[k, :], bins=bins, alpha=alph, label='$f_1$ (Sampled)', color='b')
        axis.hist(totl, bins=bins, alpha=alph, label='$f_0$ (Sampled)', color='g')
        axis.plot(meansing, powrprob * numbsamp * deltsing, label='$f_1$ (Analytic)', color='b')
        axis.plot(meandoub, convprob * numbsamp * deltdoub[0], label='$f_0$ (Numerically convolved)', color='g')
        
        axis.plot(meandoub[indxdoub], convprobpoly * numbsamp * deltdoub[indxdoub], label='$f_0$ (Fit)', color='r')
    
        axis.set_ylim([0.5, numbsamp])
        axis.set_xlabel('$f$')
        axis.set_xlim([np.amin(bins), np.amax(bins)])
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_ylabel('$N_{samp}$')
        make_legd(axis)
        plt.tight_layout()
        pathfold = os.environ["TDGU_DATA_PATH"] + '/imag/powrpdfn/'
        os.system('mkdir -p ' + pathfold)
        figr.savefig(pathfold + 'powrpdfn%04d.pdf' % n)
        plt.close(figr)
        

def savefigr(gdat, gdatmodi, figr, path):
    
    #if gdatmodi is not None and gdat.numbproc > 1:
    #    gdatmodi.lock.acquire()
    #    print 'Process %d acquiring the lock...' % gdatmodi.indxprocwork 
    
    plt.savefig(path)
    
    #if gdatmodi is not None and gdat.numbproc > 1:
    #    gdatmodi.lock.release()
    #    print 'Process %d releasing the lock...' % gdatmodi.indxprocwork 
        

def plot_elemtdim(gdat, gdatmodi, strgstat, strgmodl, strgelemtdimtype, strgelemtdimvarb, indxpoplfrst, strgfrst, \
                                                                                                          strgseco, strgtotl, strgmome='pmea', strgpdfn='post'):
    
    sizelarg = 10
    sizesmll = 1
    
    if strgstat == 'pdfn':
        legdmome = getattr(gdat, 'legd' + strgmome)
    
    limtfrst = getattr(gdat, 'limt' + strgfrst)
    limtseco = getattr(gdat, 'limt' + strgseco)
    factplotfrst = getattr(gdat, 'fact' + strgfrst + 'plot')
    factplotseco = getattr(gdat, 'fact' + strgseco + 'plot')
   
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    if strgmodl == 'fitt':
        colrtemp = gdat.fittcolrelem[indxpoplfrst]
        if strgstat == 'pdfn':
            labl = gdat.legdsampdist + ' ' + legdmome
            if strgelemtdimtype == 'bind':
                varb = getattr(gdat, strgmome + strgpdfn + strgtotl)
                varbfrst = getattr(gdat, 'bins' + strgfrst) * getattr(gdat, 'fact' + strgfrst + 'plot')
                varbseco = getattr(gdat, 'bins' + strgseco) * getattr(gdat, 'fact' + strgseco + 'plot')
                if strgtotl.startswith('hist') or strgtotl.startswith('exr') or strgtotl.startswith('incr') or np.amax(varb) <= 0.:
                    normtdim = None
                else:
                    normtdim = mpl.colors.LogNorm(0.5, vmax=np.amax(varb))
                imag = axis.pcolor(varbfrst, varbseco, varb.T, cmap='Blues', label=labl, norm=normtdim)
                make_cbar(gdat, axis, imag)
                
            else:
                if gdat.boolcondcatl:
                    varbfrst = np.zeros(gdat.numbprvlhigh)
                    varbseco = np.zeros(gdat.numbprvlhigh)
                    cntr = 0
                    for r in gdat.indxstkscond:
                        if r in gdat.indxprvlhigh:
                            varbfrst[cntr] = gdat.dictglob['poststkscond'][r][strgfrst][indxpoplfrst] * getattr(gdat, 'fact' + strgfrst + 'plot')
                            varbseco[cntr] = gdat.dictglob['poststkscond'][r][strgseco][indxpoplfrst] * getattr(gdat, 'fact' + strgseco + 'plot')
                            cntr += 1
                    axis.scatter(varbfrst, varbseco, alpha=gdat.alphelem, color=colrtemp, label=gdat.legdsamp)
        
        if strgstat == 'this' or strgstat == 'mlik':
            if strgelemtdimtype == 'bind':
                meanfrst = getattr(gdat, 'bins' + strgfrst) * getattr(gdat, 'fact' + strgfrst + 'plot')
                meanseco = getattr(gdat, 'bins' + strgseco) * getattr(gdat, 'fact' + strgseco + 'plot')
                hist = getattr(gdatmodi, strgstat + strgtotl)
                if strgtotl.startswith('hist') or strgtotl.startswith('exr') or strgtotl.startswith('incr') or np.amax(hist) <= 0.:
                    normtdim = None
                else:
                    normtdim = mpl.colors.LogNorm(0.5, vmax=np.amax(hist))
                imag = axis.pcolor(meanfrst, meanseco, hist.T, cmap='Blues', label=gdat.legdsamp, alpha=gdat.alphhist, norm=normtdim)
            else:
                varbfrst = getattr(gdatmodi, 'this' + strgfrst)[indxpoplfrst] * getattr(gdat, 'fact' + strgfrst + 'plot')
                varbseco = getattr(gdatmodi, 'this' + strgseco)[indxpoplfrst] * getattr(gdat, 'fact' + strgseco + 'plot')
                if len(varbfrst) == 0 or len(varbseco) == 0:
                    varbfrst = np.array([limtfrst[0] * factplotfrst * 0.1])
                    varbseco = np.array([limtseco[0] * factplotseco * 0.1])
                axis.scatter(varbfrst, varbseco, alpha=gdat.alphelem, color=colrtemp, label=gdat.legdsamp)
    
    # reference elements
    if strgfrst[-4:] in gdat.listnamerefr:
        strgfrsttemp = strgfrst[-4:]
    else:
        strgfrsttemp = strgfrst
    if strgseco[-4:] in gdat.listnamerefr:
        strgsecotemp = strgseco[-4:]
    else:
        strgsecotemp = strgseco
    if hasattr(gdat, 'refr' + strgfrsttemp) and hasattr(gdat, 'refr' + strgsecotemp):
        for q in gdat.indxrefr:
            if strgfrsttemp in gdat.refrliststrgfeat[q] and strgsecotemp in gdat.refrliststrgfeat[q]:
                refrvarbfrst = getattr(gdat, 'refr' + strgfrsttemp)[q] * getattr(gdat, 'fact' + strgfrsttemp + 'plot')
                refrvarbseco = getattr(gdat, 'refr' + strgsecotemp)[q] * getattr(gdat, 'fact' + strgsecotemp + 'plot')
                if len(refrvarbfrst) == 0 or len(refrvarbseco) == 0:
                    refrvarbfrst = np.array([limtfrst[0] * factplotfrst * 0.1])
                    refrvarbseco = np.array([limtseco[0] * factplotseco * 0.1])
                axis.scatter(refrvarbfrst, refrvarbseco, alpha=gdat.alphelem, color=gdat.refrcolrelem[q], label=gdat.refrlegdelem[q], s=sizelarg)

    plot_sigmcont(gdat, strgmodl, axis, strgfrst, indxpoplfrst, strgseco=strgseco)
    
    scalfrst = getattr(gdat, 'scal' + strgfrst + 'plot')
    scalseco = getattr(gdat, 'scal' + strgseco + 'plot')
    if scalfrst == 'logt':
        axis.set_xscale('log')
    if scalseco == 'logt':
        axis.set_yscale('log')
    
    axis.set_xlabel(getattr(gdat, 'labl' + strgfrst + 'totl'))
    axis.set_ylabel(getattr(gdat, 'labl' + strgseco + 'totl'))
    axis.set_xlim(np.array(limtfrst) * factplotfrst)
    axis.set_ylim(np.array(limtseco) * factplotseco)
    
    make_legd(axis)

    plt.tight_layout()
    if strgstat == 'pdfn':
        strgmometemp = strgmome
    else:
        strgmometemp = ''
    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, '%s%s' % (strgmometemp, strgtotl), nameinte=strgelemtdimvarb + 'tdim/')
    
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_sigmcont(gdat, strgmodl, axis, strgfrst, indxpoplfrst, strgseco=None):
    
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    
    if strgfrst == 'deltllik' or strgseco == 'deltllik':
        for pval in gdat.pvalcont:
            if strgfrst == 'deltllik':
                deltlliksigm = scipy.stats.chi2.ppf(1. - pval, numbcomp[indxpoplfrst])
                axis.axvline(deltlliksigm, ls='--', color='black', alpha=0.2) 
            if strgseco == 'deltllik':
                deltlliksigm = scipy.stats.chi2.ppf(1. - pval, numbcomp[indxpoplfrst])
                axis.axhline(deltlliksigm, ls='--', color='black', alpha=0.2) 
    

def plot_gene(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgydat, strgxdat, histtype='hist', \
                     indxrefrplot=None, indxydat=None, strgindxydat=None, indxxdat=None, strgindxxdat=None, plottype='none', \
                     scal=None, scalxdat=None, scalydat=None, limtxdat=None, limtydat=None, omittrue=False, nameinte='', \
                     lablxdat='', lablydat='', factxdat=1., factydat=1., histodim=False, offslegd=None, tdim=False, ydattype='totl', boolhistprio=True):
    
    if strgydat[-8:-5] == 'pop':
        boolelem = True
    else:
        boolelem = False

    if scal is None:
        if scalxdat is None:
            scalxdat = 'linr'
        if scalydat is None:
            scalydat = 'linr'
    else:
        scalxdat = scal
        scalydat = scal

    if histodim:
        figrsize = (gdat.plotsize, 0.8 * gdat.plotsize)
    else:
        figrsize = (gdat.plotsize, gdat.plotsize)

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    
    if tdim:
        xdat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgxdat, strgpdfn) * factxdat
        ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn) * factydat
    else:
        xdat = getattr(gdat, strgxdat) * factxdat
        if histtype == 'histptfn':
            ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'histptfn' + strgydat[4:], strgpdfn) * factydat
        else:
            ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn) * factydat
    
    if indxxdat is not None:
        xdat = xdat[indxxdat]
    if indxydat is not None:
        ydat = ydat[indxydat]
    
    xerr = np.zeros((2, xdat.size))
    
    if tdim:
        axis.scatter(xdat, ydat, alpha=gdat.alphelem, color=colr, label=gdat.legdsamp)
    else:
        if histodim:
            # temp
            if strgxdat[4:] in gdat.fittliststrgfeattotl:
                deltxdat = getattr(gdat, 'delt' + strgxdat[4:]) * factxdat
                binsxdat = getattr(gdat, 'bins' + strgxdat[4:]) * factxdat
            else:
                deltxdat = getattr(gdat, 'delt' + strgxdat[4:]) * factxdat
                binsxdat = getattr(gdat, 'bins' + strgxdat[4:]) * factxdat

            xdattemp = binsxdat[:-1] + deltxdat / 2.
   
    if strgmodl == 'fitt':
        if boolelem:
            if strgydat.startswith('cmpl'):
                legd = gdat.fittlegdelem[int(strgydat[-5])]
                colr = gdat.fittcolrelem[int(strgydat[-5])]
            else:
                legd = gdat.fittlegdelem[int(strgydat[-1])]
                colr = gdat.fittcolrelem[int(strgydat[-1])]
        else:
            legd = gdat.fittlegd
            colr = gdat.fittcolr
        
        if strgstat == 'pdfn':
            if histtype == 'histptfn':
                yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'histptfn' + strgydat[4:], strgpdfn, strgmome='errr') * factydat
            else:
                yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn, strgmome='errr') * factydat
            if indxydat is not None:
                yerr = yerr[[slice(None)] + indxydat]
            
            # label
            if strgydat.startswith('hist'):
                ##  element distribution
                labl = gdat.legdsampdist
            else:
                ##  other
                labl = gdat.legdsampdist
            
            # draw points
            indxerrr = np.where((yerr[0, :] > 0.) | (yerr[1, :] > 0.))[0]
            if indxerrr.size > 0:
                labltemp = None
            else:
                labltemp = labl
            temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, xerr=xerr, label=labl, \
                                                                                                 marker='o', ls='', markersize=5, color=colr, lw=1, capsize=5)

            # draw error-bar caps 
            if indxerrr.size > 0:
                temp, listcaps, temp = axis.errorbar(xdat[indxerrr], ydat[indxerrr], yerr=yerr[:, indxerrr], xerr=xerr[:, indxerrr], \
                                                                                                marker='o', ls='', markersize=5, color=colr, lw=1, capsize=5)
                for caps in listcaps:
                    caps.set_markeredgewidth(1)

        elif strgstat == 'this' or strgstat == 'mlik':
            
            if strgstat == 'this':
                legd = gdat.legdsamp
            else:
                legd = gdat.legdmlik

            if histodim:
                axis.bar(xdattemp, ydat, deltxdat, label=gdat.legdsamp, alpha=0.5, linewidth=1, edgecolor=colr)
            else:
                if plottype == 'errr':
                    yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgydat, strgpdfn, strgmome='errr') * factydat

                    if indxydat is not None:
                        yerr = yerr[[slice(None)] + indxydat]
                    temp, listcaps, temp = axis.errorbar(xdat, ydat, yerr=yerr, xerr=xerr, \
                                            marker='o', ls='', markersize=5, label=legd, lw=1, capsize=5, color=colr)
                    for caps in listcaps:
                        caps.set_markeredgewidth(1)
                else:
                    axis.plot(xdat, ydat, label=gdat.legdsamp, alpha=0.5, color=colr)
    
    # reference histogram
    if not omittrue:
        for q in gdat.indxrefr:
            
            if strgydat.startswith('psfn'):
                name = 'psfnexpr'
            else:
                if boolelem:
                    if strgydat[-12:-8] in gdat.listnamerefr:
                        name = 'refr' + strgydat[:-12] + 'pop%d' % q + strgydat[-4:]
                    else:
                        name = 'refr' + strgydat[:-8] + 'pop%d' % q + strgydat[-4:]
                else:
                    name = 'refr' + strgydat
            
            if not hasattr(gdat, name):
                continue

            ydattemp = getattr(gdat, name)
            
            ydat = ydattemp * factydat
            if indxydat is not None:
                ydat = ydat[indxydat]
            
            if strgydat[-8:-5] == 'pop':
                legd = gdat.refrlegdelem[q]
                colr = gdat.refrcolrelem[q]
            else:
                legd = gdat.refrlegd
                colr = gdat.refrcolr
    
            if histodim:
                axis.bar(xdattemp, ydat, deltxdat, color=colr, label=legd, alpha=gdat.alphhist, linewidth=1, edgecolor=colr)
            else:
                axis.plot(xdat, ydat, color=colr, label=legd, alpha=gdat.alphline)
                           
            try:
                if histodim:
                    if histtype == 'histptfn':
                        ptfn = getattr(gdat, 'trueptfn' + strgydat[4:])
                    axis.plot(xdattemp, 10. * ptfn, color='purple', label='PTFN', alpha=gdat.alphline)
            except:
                pass

            if not boolelem:
                break
    
    # external reference histogram
    if histodim and strgydat == 'histfluxpop0':
        try:
            if gdat.listprefhistfluxlabl is not None:
                for k in range(len(gdat.listprefhistfluxlabl)):
                    if gdat.listprefhistfluxtype[k] == 'shad':
                        axis.plot(gdat.listprefhistfluxflux[k][0], gdat.listprefhistfluxhist[k][0], color='m', label=gdat.listprefhistfluxlabl[k])
                        enerpoly = np.empty(gdat.listprefhistfluxflux[k][1].size + gdat.listprefhistfluxflux[k][2].size)
                        enerpoly[:gdat.listprefhistfluxflux[k][1].size] = gdat.listprefhistfluxflux[k][1]
                        enerpoly[gdat.listprefhistfluxflux[k][1].size:] = gdat.listprefhistfluxflux[k][2][::-1]
                        sbrtpoly = np.empty(gdat.listprefhistfluxflux[k][1].size + gdat.listprefhistfluxflux[k][2].size)
                        sbrtpoly[:gdat.listprefhistfluxflux[k][1].size] = gdat.listprefhistfluxhist[k][1]
                        sbrtpoly[gdat.listprefhistfluxflux[k][1].size:] = gdat.listprefhistfluxhist[k][2][::-1]
                        axis.fill(enerpoly, sbrtpoly, color='m', alpha=0.5)
                    else:
                        axis.errorbar(gdat.listprefhistfluxflux[k], gdat.listprefhistfluxhist[k], label=gdat.listprefhistfluxlabl[k], color='m')
        except:
            pass

    if strgydat.startswith('histcntp'):
        if gdat.numbpixl > 1:
            ydattemp = getattr(gdat, 'histcntpdata' + strgydat[-8:])
        else:
            ydattemp = getattr(gdat, 'histcntpdata' + strgydat[-4:])
        axis.bar(xdattemp, ydattemp, deltxdat, color='black', label='Data', alpha=gdat.alphhist, linewidth=1, edgecolor='black')
                
    # axis scales
    if scalxdat == 'logt':
        axis.set_xscale('log')
    if scalydat == 'logt':
        if np.where(ydat > 0.)[0].size > 0:
            axis.set_yscale('log')
    
    # axis labels
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)

    # superimpose prior on the feature
    liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
    liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
    
    ptch = None
    line = None

    if strgydat.startswith('hist') and strgydat != 'histdefl' and strgydat != 'histdeflelem' and boolhistprio:
        if strgydat[-8:-5] == 'pop':
            strgtemp = strgydat[4:-8]
            if strgtemp in liststrgfeatprio[int(strgydat[-5])]:
                xdatprio = getattr(gdat, strgmodl + strgxdat + 'prio') * factxdat
                if gdat.datatype == 'mock' and not omittrue:
                    for q in gdat.indxrefr:
                        if gdat.refrnumbelempopl[q] == 0:
                            continue
                        if strgtemp in gdat.trueliststrgfeatprio[q]:
                            truexdatprio = getattr(gdat, 'true' + strgxdat + 'prio') * factxdat
                            trueydatsupr = getattr(gdat, 'true' + strgydat + 'prio') * factydat
                            trueydatsupr = retr_fromgdat(gdat, gdatmodi, strgstat, 'true', strgydat + 'prio', strgpdfn) * factydat
                            axis.plot(truexdatprio, trueydatsupr, ls='-', alpha=gdat.alphline, color=gdat.refrcolrelem[q])

                if strgmodl != 'true':
                    ydatsupr = retr_fromgdat(gdat, gdatmodi, strgstat, 'fitt', strgydat + 'prio', strgpdfn) * factydat
                    if strgstat == 'pdfn':
                        yerrsupr = retr_fromgdat(gdat, gdatmodi, strgstat, 'fitt', strgydat + 'prio', strgpdfn, strgmome='errr') * factydat
                        labl = gdat.legdsampdist + ' hyper-distribution'
                        ptch, line = tdpy.util.plot_braz(axis, xdatprio, ydatsupr, yerr=yerrsupr, lcol='lightgrey', dcol='grey', labl=labl)
                    else:
                        axis.plot(xdatprio, ydatsupr, ls='--', alpha=gdat.alphline, color=gdat.fittcolrelem[int(strgydat[-5])])
   
    for name, valu in gdat.__dict__.iteritems():
        if name.startswith('refrplot'):
            if name[8:12] == 'hist' and name[12:16] == strgydat[4:] and name[16:19] == 'pop' and int(name[-1]) == indxpopltemp:
                colr = getattr(gdat, name + 'colr')
                linestyl = getattr(gdat, name + 'linestyl')
                axis.plot(valu[0, :], valu[1, :], ls=linestyl, color=colr)

    if strgydat.startswith('hist') and strgydat[4:-8] == 'deltllik':
        plot_sigmcont(gdat, strgmodl, axis, strgxdat[4:], int(strgydat[-1]))
   
    if indxydat is not None:
        strgydat += strgindxydat
    
    if indxxdat is not None:
        strgxdat += strgindxxdat
    
    if limtxdat is not None:
        axis.set_xlim(limtxdat)
    else:
        axis.set_xlim([np.amin(xdat), np.amax(xdat)])
    if limtydat is not None:
        axis.set_ylim([limtydat[0] * factydat, limtydat[1] * factydat])
    else:
        axis.set_ylim([np.amin(ydat), np.amax(ydat)])
    
    if ydattype != 'totl':
        strgydat += ydattype

    try:
        make_legd(axis, offs=offslegd, ptch=ptch, line=line)
    except:
        print 'Legend failed when'
        print 'strgstat'
        print strgstat
        print 'strgmodl'
        print strgmodl
        print 'strgydat'
        print strgydat
        print

    plt.tight_layout()
    if histtype == 'histptfn':
        path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, 'histptfn' + strgydat[4:], nameinte=nameinte)
    else:
        path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, strgydat, nameinte=nameinte)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def plot_scatassc(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, q, l, strgfeat, plotdiff=False):
    
    if plotdiff:
        figrsize = (gdat.plotsize, 0.7 * gdat.plotsize)
    else:
        figrsize = (gdat.plotsize, gdat.plotsize)
    figr, axis = plt.subplots(1, 1, figsize=figrsize)
    
    # prepare data to be plotted
    xdat = np.copy(getattr(gdat, 'refr' + strgfeat)[q][0, :])
    xerr = tdpy.util.retr_errrvarb(getattr(gdat, 'refr' + strgfeat)[q])
   
    minmplot = getattr(gdat, 'minm' + strgfeat)
    maxmplot = getattr(gdat, 'maxm' + strgfeat)
    binsplot = getattr(gdat, 'bins' + strgfeat)
    factplot = getattr(gdat, 'fact' + strgfeat + 'plot')

    ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgfeat + 'asscpop%dpop%d' % (q, l), strgpdfn) * factplot
    
    yerr = np.zeros((2, ydat.size))
    if strgstat == 'pdfn':
        yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgfeat + 'asscpop%dpop%d' % (q, l), strgpdfn, strgmome='errr') * factplot
    
    xdat *= factplot
    xerr *= factplot
    if plotdiff:
        ydat = 100. * (ydat - xdat) / xdat
    
    # handle the case when there is a single reference element
    if yerr.ndim == 1:
        ydat = np.array([ydat])
        yerr = yerr[:, None]
    
    # plot all associations
    if plotdiff:
        indx = np.where(ydat > -100.)[0]
    else:
        indx = np.where(ydat > 0.)[0]
    if indx.size > 0:
        axis.errorbar(xdat[indx], ydat[indx], ls='', yerr=yerr[:, indx], xerr=xerr[:, indx], lw=1, marker='o', markersize=5, color='black')
    
    # temp -- plot associations inside the comparison area
    if plotdiff:
        axis.axhline(0., ls='--', alpha=gdat.alphline, color='black')
    else:
        axis.plot(factplot * binsplot, factplot * binsplot, ls='--', alpha=gdat.alphline, color='black')
    
    lablxdat = getattr(gdat, 'labl' + strgfeat + 'refr')
    lablydat = getattr(gdat, 'labl' + strgfeat + 'samp')
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)
    boollogtxaxi = False
    boollogtyaxi = False
    if indx.size > 0 and getattr(gdat, 'scal' + strgfeat + 'plot') == 'logt':
        if not plotdiff:
            axis.set_yscale('log')
            boollogtyaxi = True
        axis.set_xscale('log')
        boollogtaxis = True
   
    if plotdiff:
        limtydat = np.array([-100., 100.])
    else:
        limtydat = factplot * np.array([minmplot, maxmplot])
    limtxdat = [factplot * minmplot, factplot * maxmplot]
    
    # overplot text
    if 'etag' in gdat.refrliststrgfeat[q]:
        for k in range(indx.size):
            if boollogtxaxi:
                sizexoff = 0.01 * xdat[indx[k]]
            else:
                sizexoff = 0.01 * (limtxdat[1] - limtxdat[0])
            if boollogtyaxi:
                sizeyoff = 0.01 * ydat[indx[k]]
            else:
                sizeyoff = 0.01 * (limtydat[1] - limtydat[0])
            axis.text(xdat[indx[k]] + sizexoff, ydat[indx[k]] + sizeyoff, gdat.refretag[q][indx[k]], verticalalignment='center', horizontalalignment='center', \
                                                                                                                                                        color='red', fontsize=1)

    axis.set_ylim(limtydat)
    axis.set_xlim(limtxdat)
   
    plt.tight_layout()
    if plotdiff:
        strgtype = 'diff'
    else:
        strgtype = ''
    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, 'scatassc' + strgfeat + '%spop%dpop%d' % (strgtype, q, l), nameinte='assc/')
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)


def plot_scatcntp(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, indxevttplot, indxenerplot=None):
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    ydat = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'cntpmodl', strgpdfn)
    if indxenerplot is None:
        xdat = gdat.cntpdata[:, :, indxevttplot].flatten()
        ydat = ydat[:, :, indxevttplot].flatten()
        nameplot = 'scatcntpevt%d' % (indxevttplot)
        if strgstat == 'pdfn':
            indxvarb = [slice(None), slice(None), indxevttplot]
    else:
        xdat = gdat.cntpdata[indxenerplot, :, indxevttplot]
        ydat = ydat[indxenerplot, :, indxevttplot]
        nameplot = 'scatcntpen%02devt%d' % (indxenerplot, indxevttplot)
        if strgstat == 'pdfn':
            indxvarb = [indxenerplot, slice(None), indxevttplot]
    if strgstat == 'pdfn':
        yerr = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, 'cntpmodl', strgpdfn, strgmome='errr', indxvarb=indxvarb)
    colr = gdat.fittcolr

    if strgstat == 'pdfn':
        axis.errorbar(xdat, ydat, yerr=yerr, marker='o', ls='', markersize=5, color=gdat.fittcolr, capsize=5)
    else:
        axis.plot(xdat, ydat, marker='o', ls='', markersize=5, color=gdat.fittcolr)
    gdat.limtcntpdata = [gdat.binscntpdata[0], gdat.binscntpdata[-1]]
    axis.set_xlim(gdat.limtcntpdata)
    axis.set_ylim(gdat.limtcntpdata)
    axis.set_ylabel('$k^{modl}$')
    axis.set_xlabel('$k^{data}$')
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.tight_layout()

    path = retr_plotpath(gdat, gdatmodi, strgpdfn, strgstat, strgmodl, nameplot)
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    
    
def plot_indxprox(gdat):

    numbbins = 40
    numbfluxprox = len(gdat.indxpixlprox)
    bins = np.empty((gdat.numbprox, numbbins + 1))
    indxpixlproxsize = np.empty((numbfluxprox, gdat.numbpixlfull))
    for h in range(gdat.numbprox):
        for j in gdat.indxpixlfull:
            try:
                indxpixlproxsize[h, j] = gdat.indxpixlprox[h][j].size
            except:
                indxpixlproxsize[h, j] = gdat.numbpixlfull
        bins[h, :] = np.logspace(np.log10(np.amin(indxpixlproxsize[h, :])), np.log10(np.amax(indxpixlproxsize[h, :])), numbbins + 1)
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for h in range(gdat.numbprox):
        axis.hist(indxpixlproxsize[h, :], bins=bins[h, :], log=True, label='Flux bin %d' % h, alpha=gdat.alphhist)
    axis.set_xscale('log')
    axis.axvline(gdat.numbpixlfull, label='ROI', ls='--')
    axis.set_xlabel('Number of pixels')
    axis.set_ylabel("Number of tables")
    make_legd(axis)
    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'init/indxprox.pdf')
    plt.close()
    
    
def plot_psfn_type():
    
    devi = np.linspace(0., 5., 100)
    y = np.zeros((x.size, 5))

    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
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

    make_legd(axis)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylim([1e-3, None])
    plt.show()
    
    
def plot_evidtest():
    
    minmgain = -1.
    maxmgain = 5.
    minmdevi = 0.
    maxmdevi = 5.
    gain = np.linspace(minmgain, maxmgain, 100)
    devi = np.linspace(minmdevi, maxmdevi, 100)

    evid = np.log(np.sqrt(1. + np.exp(2. * gain[None, :])) * np.exp(-devi[:, None]**2 / 2. / (1. + 1. / np.exp(2. * gain[None, :]))))
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    figr.suptitle('Log-Bayesian Evidence For Lower-Dimension Model', fontsize=18)
    imag = axis.imshow(evid, extent=[minmgain, maxmgain, minmdevi, maxmdevi], cmap='winter', origin='lower')
    cset1 = plt.contourf(gain, devi, evid, cmap='winter')
    axis.set_xlabel('Information gain')
    axis.set_ylabel('Goodness of fit')
    plt.colorbar(imag, ax=axis, fraction=0.03)

    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'evidtest.pdf')
    plt.close(figr)
    
    
def plot_histlgalbgalelemstkd(gdat, strgpdfn, indxpoplplot, strgbins, strgfeat=None):
    
    if strgfeat is not None:
        numbparaplot = gdat.numbbinsplot
    else:
        numbparaplot = 1

    if strgbins == 'cumu':
        numbrows = 1
        numbcols = 1
    else:
        numbcols = 2
        if strgbins == 'full':
            numbrows = numbparaplot / 2
        else:
            numbrows = 2
    
    if strgfeat is not None:
        indxfeatsign = gdat.fittliststrgfeatsign[indxpoplplot].index(strgfeat)
    else:
        indxfeatsign = np.arange(len(gdat.fittliststrgfeatsign))
    
    histlgalbgalelemstkd = getattr(gdat, strgpdfn + 'histlgalbgalelemstkd')

    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize), sharex='all', sharey='all')
    if numbrows == 1:
        axgr = [axgr]            
    for a, axrw in enumerate(axgr):
        if numbcols == 1:
            axrw = [axrw]
        for b, axis in enumerate(axrw):
            if strgfeat is not None:
                h = a * 2 + b
                if strgbins == 'full':
                    indxlowr = h
                    indxuppr = h + 1
                elif strgbins == 'cumu':
                    indxlowr = 0
                    indxuppr = numbparaplot
                else:
                    if h < 3:
                        indxlowr = 2 * h
                        indxuppr = 2 * (h + 1)
                    else:
                        indxlowr = 2 * h
                        indxuppr = numbparaplot
                temp = np.sum(histlgalbgalelemstkd[indxpoplplot][:, :, indxlowr:indxuppr, indxfeatsign], 2).T
            else:
                temp = np.sum(np.sum(histlgalbgalelemstkd[indxpoplplot], 2), 2).T
                
            if np.where(temp > 0.)[0].size > 0:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi, norm=mpl.colors.LogNorm(vmin=0.5, vmax=None))
            else:
                imag = axis.imshow(temp, interpolation='nearest', origin='lower', cmap='BuPu', extent=gdat.exttrofi)
                
            if strgfeat is not None:
                bins = getattr(gdat, 'bins' + strgfeat)
            
            # superimpose reference elements
            for q in gdat.indxrefr:
                if gdat.refrnumbelem[q] == 0:
                    continue
                # temp -- backcomp
                try:
                    reframpl = getattr(gdat, 'refr' + gdat.refrlistnamefeatampl[q])
                except:
                    reframpl = getattr(gdat, 'refr' + gdat.listnamefeatamplrefr[q])
                if strgfeat in gdat.refrliststrgfeat[q]:
                    refrfeat = getattr(gdat, 'refr' + strgfeat)[q]
                    if len(refrfeat) > 0:
                        indxelem = np.where((bins[indxlowr] < refrfeat[0, :]) & (refrfeat[0, :] < bins[indxuppr]))[0]
                    else:
                        indxelem = np.array([])
                else:
                    indxelem = np.arange(gdat.refrnumbelem[q])
                # temp -- backcomp
                try:
                    mrkrsize = retr_mrkrsize(gdat, reframpl[q][0, indxelem], gdat.refrlistnamefeatampl[q])
                except:
                    mrkrsize = retr_mrkrsize(gdat, reframpl[q][0, indxelem], gdat.listnamefeatamplrefr[q])

                if indxelem.size > 0:
                    axis.scatter(gdat.anglfact * gdat.refrlgal[q][0, indxelem], gdat.anglfact * gdat.refrbgal[q][0, indxelem], \
                                                s=mrkrsize, alpha=gdat.alphelem, marker=gdat.refrlistelemmrkrhits[q], lw=2, color=gdat.refrcolrelem[q])

            if a == numbrows - 1:
                axis.set_xlabel(gdat.labllgaltotl)
            else:
                axis.set_xticklabels([])
            if b == 0:
                axis.set_ylabel(gdat.lablbgaltotl)
            else:
                axis.set_yticklabels([])

            draw_frambndr(gdat, axis)
            
            if strgbins != 'cumu':
                lablfeat = getattr(gdat, 'labl' + strgfeat)
                factfeatplot = getattr(gdat, 'fact' + strgfeat + 'plot')
                titl = tdpy.util.mnp.exp(factfeatplot * bins[indxlowr]) + ' < $%s$ < ' % lablfeat + tdpy.util.mnp.exp(factfeatplot * bins[indxuppr])
                axis.set_title(titl)
    
    if strgfeat is not None:
        lablfeattotl = getattr(gdat, 'labl' + strgfeat + 'totl')
        plt.figtext(0.5, 0.95, '%s' % lablfeattotl, ha='center', va='center')
    axiscomm = figr.add_axes([0.87, 0.2, 0.02, 0.6])
    cbar = figr.colorbar(imag, cax=axiscomm)

    plt.subplots_adjust()
    #plt.subplots_adjust(left=0.18, top=.9, right=0.82, bottom=0.15, hspace=0.08, wspace=0.08)
    if strgbins == 'cumu':
        strgtemp = ''
    else:
        strgtemp = strgfeat
    path = getattr(gdat, 'path' + strgpdfn + 'finl') + 'histlgalbgalelemstkd%s%spop%d' % (strgbins, strgtemp, indxpoplplot) + '.pdf'
    figr.savefig(path)
    plt.close(figr)
       
    
def plot_king(gdat):

    angl = rad2deg(gdat.binsangl)

    figr, axgr = plt.subplots(1, 2, figsize=(2 * gdat.plotsize, gdat.plotsize))
    figr.suptitle('King Function', fontsize=20)
    for k, axis in enumerate(axgr):
        if k == 0:
            sigmlist = [0.25]
            gammlist = [1.01, 2.5, 10.]
        else:
            sigmlist = [0.1, 0.25, 1.]
            gammlist = [2.]
        for sigm in sigmlist:
            for gamm in gammlist:
                axis.plot(angl, retr_singking(angl, sigm, gamm), label=r'$\sigma = %.4g, \gamma = %.3g$' % (sigm, gamm))
        make_legd(axis)
        axis.set_yscale('log')
        axis.set_xlabel(gdat.lablgangtotl)
        axis.set_xlabel(r'$\mathcal{K}$')
        
    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'king.pdf')
    plt.close(figr)
    
   
def plot_intr(gdat):
    
    if gdat.verbtype > 0:
        print 'Making PCAT introductory plots...'

    #plot_grap(plottype='meta', verbtype=1)
    plot_grap(plottype='lght0000', verbtype=1)
    #plot_grap(plottype='lght0001', verbtype=1)
    #plot_grap(plottype='lght0002', verbtype=1)
    #plot_grap(plottype='lght0003', verbtype=1)
    #plot_grap(plottype='lens0000', verbtype=1)
    plot_grap(plottype='lens0001', verbtype=1)
    
    with plt.xkcd():

        from matplotlib import patheffects
        mpl.rcParams['path.effects'] = [patheffects.withStroke(linewidth=0)]

        figr, axis = plt.subplots(figsize=(2 * gdat.plotsize, gdat.plotsize))

        catl = np.arange(80)
        probcatl = pss.pmf(catl, 30.) + 0.5 * pss.pmf(catl, 60.)
        axis.plot(catl, probcatl)
        axis.set_xticks([10, 30, 60])
        axis.set_xticklabels(["Crackpot's Catalog", "Best-fit catalog", "Not-so-best-fit catalog"])
        axis.set_yticks([])
        text = axis.set_title("Exploring the catalog space with Probabilistic cataloging")
        text.set_position([.5, 1.05])
        axis.set_xlabel('Catalog index')
        axis.set_ylabel("Probability")
        
        axis.tick_params(axis='x', colors='#B6E954')
        axis.tick_params(axis='y', colors='#B6E954')
        axis.spines['bottom'].set_color('#B6E954')
        axis.spines['top'].set_color('#B6E954') 
        axis.spines['right'].set_color('#B6E954')
        axis.spines['left'].set_color('#B6E954')
        axis.yaxis.label.set_color('#B6E954')
        axis.xaxis.label.set_color('#B6E954')
        axis.title.set_color('#B6E954')

        axis.set_axis_bgcolor('black')
        figr.set_facecolor('black')
        plt.tight_layout()
        figr.savefig(gdat.pathimag + 'talkintr.pdf', facecolor=figr.get_facecolor())
        plt.close()  
        
        
def plot_eval(gdat, indxpoplplot):

    psfntemp = np.copy(gdat.psfnexpr[0, :, 0])
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    for k in range(gdat.numbprox + 1):
        if k == 0 or k == gdat.numbprox:
            alph = 1.
            colr = gdat.fittcolrelem[indxpoplplot]
            if k == 0:
                labl = 'Dimmest PS'
            else:
                labl = 'Brightest PS'
        else:
            alph = 0.2
            labl = None
            colr = 'black'
        axis.plot(gdat.binsangl * gdat.anglfact, gdat.binsprox[k] * psfntemp, label=labl, color=colr, alpha=alph)
        axis.set_xlim([np.amin(gdat.binsangl) * gdat.anglfact, np.amax(gdat.binsangl) * gdat.anglfact])
        if k > 0:
            axis.axvline(gdat.anglfact * gdat.maxmangleval[k-1], ls='--', alpha=alph, color=colr)
    axis.set_yscale('log')
    axis.set_xlabel(gdat.lablgangtotl)
    axis.set_ylabel(gdat.lablsbrttotl)
    
    limt = gdat.specfraceval * np.amax(gdat.binsprox[0] * psfntemp)

    if limt != 0.:
        axis.axhline(limt, color='red', ls=':', label='Flux floor')
    
    #axis.set_ylim([1e-3 * limt, None])
    #maxmangltemp = np.interp(1e-1 * limt, gdat.binsprox[k] * psfntemp[::-1], gdat.binsangl[::-1] * gdat.anglfact)
    #axis.set_xlim([None, maxmangltemp])
    
    make_legd(axis)

    plt.tight_layout()
    figr.savefig(gdat.pathinit + 'eval.pdf')
    plt.close(figr)


def plot_mosa(gdat, strgpdfn):

    # empty global object
    gdatmodi = tdpy.util.gdatstrt()
    
    listsamp = getattr(gdat, 'list' + strgpdfn + 'samp')
    listsampunit = getattr(gdat, 'list' + strgpdfn + 'sampunit')

    numbrows = 3
    numbcols = 2
    numbsampmosa = numbrows * numbcols
    if numbsampmosa <= gdat.numbsamptotl:
        indxsampmosa = choice(gdat.indxsamptotl, size=numbsampmosa, replace=False)
        for l in gdat.fittindxpopl:
            for i in gdat.indxener:
                for m in gdat.indxevttplot:
                    
                    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
                    for a, axrw in enumerate(axgr):
                        for b, axis in enumerate(axrw):
                            
                            n = indxsampmosa[numbcols*a+b]
                            gdatmodi.thissamp = listsamp[n, :].flatten()
                            gdatmodi.thissampunit = listsampunit[n, :].flatten()
                            if gdat.fittnumbtrap > 0:
                                gdatmodi.thisindxelemfull = getattr(gdat, 'list' + strgpdfn + 'indxelemfull')[n]
                                proc_samp(gdat, gdatmodi, 'this', 'fitt', boolinit=True)

                            if a == numbrows - 1:
                                axis.set_xlabel(gdat.labllgaltotl)
                            else:
                                axis.set_xticklabels([])
                            if b == 0:
                                axis.set_ylabel(gdat.lablbgaltotl)
                            else:
                                axis.set_yticklabels([])
                            
                            imag = retr_imag(gdat, axis, gdat.cntpdata, '', 'fitt', 'cntpdata', i, m)
                            supr_fram(gdat, gdatmodi, 'this', 'fitt', axis, l)
                    
                    if gdat.enerbins:
                        plt.figtext(0.5, 0.93, gdat.strgener[i], ha='center', va='center')
                    axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
                    cbar = figr.colorbar(imag, cax=axiscomm)
                    cbar.set_ticks(gdat.tickcntpdata)
                    cbar.set_ticklabels(gdat.lablcntpdata)
                    plt.subplots_adjust()
                    #plt.subplots_adjust(left=0.1, top=.91, hspace=0.03, wspace=0.1, bottom=0.09)
                    if l == 1:
                        strg = ''
                    else:
                        strg = 'pop%d' % l
                    pathfinl = getattr(gdat, 'path' + strgpdfn + 'finl')
                    if m is None:
                        path = pathfinl + 'mosa' + strg + 'en%02dA.pdf' % (gdat.indxenerincl[i])
                    else:
                        path = pathfinl + 'mosa' + strg + 'en%02devtt%d.pdf' % (gdat.indxenerincl[i], gdat.indxevttincl[m])
                    figr.savefig(path)
                    plt.close(figr)
    else:
        if gdat.verbtype > 0:
            print 'Skipping the mosaic plot...'


def plot_grap(plottype, verbtype=0):
        
    figr, axis = plt.subplots(figsize=(6, 6))

    grap = nx.DiGraph()
    if plottype == 'meta':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'magenta']


    if plottype == 'lens0001':
        listcolr = ['olive', 'olive', 'black', 'magenta', 'magenta', 'magenta', 'magenta', 'magenta', 'olive', 'olive', 'olive', 'olive', 'olive', \
                                                                                                                                        r'black', 'olive', 'black']

    if plottype == 'lght0000':
        listcolr = [r'olive', r'black', r'magenta', r'magenta', 'magenta', r'magenta', r'olive', r'olive', r'black', r'olive', r'olive', r'black', r'olive']
    



    if plottype == 'lght0001':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta', 'black']

    if plottype == 'lght0002':
        listcolr = ['black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', 'olive', 'olive', 'olive', 'magenta', \
                                                                                                    'magenta', 'magenta', 'magenta', 'magenta', 'black']
    if plottype == 'lght0003':
        listcolr = ['black', 'black', 'black', 'olive', 'black', 'olive', 'olive', 'black', 'olive', \
                                                                                                    'olive', 'olive', 'magenta', 'magenta', 'magenta', 'magenta']
    
    if plottype == 'lens0000':
        listcolr = ['olive', 'black', 'black', 'olive', 'olive', 'olive', 'olive', 'black', 'olive', 'magenta', 'magenta', 'magenta']


    if plottype.startswith('meta'):
        grap.add_edges_from([ \
                             ('meanelem', 'numbelem'), \
                             ('modl','data'), \
                             ('psfp', 'modl'), \
                             ('feat','modl'), \
                             ('numbelem','feat'), \
                             ('ampldistslop', 'ampl'), \
                            ])
    
    if plottype.startswith('lght') or plottype.startswith('lens'):
        grap.add_edges_from([ \
                             ('meanelem', 'numbelem'), \
                             ('modl','data'), \
                             ('psfp', 'modl'), \
                             ('bacp', 'modl'), \
                             ('lgal','modl'), \
                             ('bgal','modl'), \
                             ('numbelem','lgal'), \
                             ('numbelem','bgal'), \
                            ])
    
    if plottype.startswith('lght'):
        grap.add_edges_from([ \
                             ('ampldistslop', 'ampl'), \
                             ('ampl', 'modl'), \
                             ('numbelem','ampl'), \
                             ('numbelem', 'sind'), \
                             ('sind','modl'), \
                            ])
    
    if plottype.startswith('lens'):
        grap.add_edges_from([ \
                             ('lenp', 'modl'), \
                             ('defsdistslop', 'defs'), \
                             ('defs', 'modl'), \
                             ('numbelem','defs'), \
                            ])
    
    if plottype == 'lens0001':
        grap.add_edges_from([ \
                             ('asca', 'modl'), \
                             ('numbelem','asca'), \
                             ('acut', 'modl'), \
                             ('numbelem','acut'), \
                            ])
    
    if plottype == 'lght0001' or plottype == 'lght0002':
        grap.add_edges_from([ \
                             ('sinddistmean', 'sind'), \
                            ])
    
    if plottype == 'lght0002':
        grap.add_edges_from([ \
                             ('numbelem', 'expc'), \
                             ('expc', 'modl'), \
                            ])
    
    if plottype == 'lght0003':
        grap.add_edges_from([ \
                             ('spatdistcons', 'lgal'), \
                             ('spatdistcons', 'bgal'), \
                            ])
        
    labl = {}
    if plottype.startswith('lens'):
        nameelem = r'\rm{sub}'
    else:
        nameelem = r'\rm{pts}'
    if plottype.startswith('lght') and (plottype == 'lght0001' or plottype == 'lght0002'):
        labl['numbelem'] = r'$\vec{N}_{%s}$' % nameelem
        labl['meanelem'] = r'$\vec{\mu}_{%s}$' % nameelem
    else:
        labl['numbelem'] = '$N_{%s}$' % nameelem
        labl['meanelem'] = r'$\mu_{%s}$' % nameelem
    
    if plottype.startswith('lght'):
        if plottype == 'lght0000' or plottype == 'lght0003':
            labl['ampldistslop'] = r'$\alpha$'
        else:
            labl['ampldistslop'] = r'$\vec{\alpha}$'
    if plottype.startswith('lens'):
        labl['defsdistslop'] = r'$\beta$'
    
    if plottype == 'lght0001' or plottype == 'lght0002':
        labl['sinddistmean'] = r'$\vec{\beta}$'
    
    if plottype == 'lght0003':
        labl['spatdistcons'] = r'$\gamma$'
    if plottype.startswith('lens'):
        labl['lenp'] = r'$\vec{\chi}$'
    labl['psfp'] = r'$\vec{\eta}$'
    labl['bacp'] = r'$\vec{A}$'
    labl['lgal'] = r'$\vec{\theta_1}$'
    labl['bgal'] = r'$\vec{\theta_2}$'
    if plottype.startswith('meta'):
        labl['feat'] = r'$\vec{\xi}$'
    else:
        if plottype.startswith('lght'):
            labl['sind'] = r'$\vec{s}$'
            labl['ampl'] = r'$\vec{f}$'
        else:
            labl['defs'] = r'$\vec{\alpha_{\rm{s}}}$'
    if plottype == 'lens0001':
        labl['asca'] = r'$\vec{\theta_{\rm{s}}}$'
        labl['acut'] = r'$\vec{\theta_{\rm{c}}}$'
        
    if plottype == 'lght0002':
        labl['expc'] = r'$\vec{E_{\rm{c}}}$'
    labl['modl'] = r'$M_D$'
    labl['data'] = r'$D$'
    
    posi = nx.circular_layout(grap)
    posi['sinddistmean'] = np.array([0.4, 0.15])
    if plottype == 'lght0003':
        posi['spatdistcons'] = np.array([-0.2, 0.15])
    if plottype.startswith('lght'):
        posi['numbelem'] = np.array([0., 0.075])
        posi['meanelem'] = np.array([0., 0.15])
        posi['ampldistslop'] = np.array([0.2, 0.15])
    if plottype.startswith('lens'):
        posi['numbelem'] = np.array([-0.1, 0.075])
        posi['meanelem'] = np.array([-0.1, 0.15])
        posi['defsdistslop'] = np.array([0.1, 0.15])
    
    if plottype.startswith('lght'):
        if plottype == 'lght0002':
            posi['psfp'] = np.array([0.7, -0.0])
            posi['bacp'] = np.array([0.9, -0.0])
        else:
            posi['psfp'] = np.array([0.5, -0.0])
            posi['bacp'] = np.array([0.7, -0.0])
    if plottype == 'lens0000':
        posi['psfp'] = np.array([0.3, -0.0])
        posi['bacp'] = np.array([0.5, -0.0])
        posi['lenp'] = np.array([0.7, -0.0])
    if plottype == 'lens0001':
        posi['psfp'] = np.array([0.7, -0.0])
        posi['bacp'] = np.array([0.9, -0.0])
        posi['lenp'] = np.array([1.1, -0.0])
    posi['lgal'] = np.array([-0.3, -0.0])
    posi['bgal'] = np.array([-0.1, -0.0])
    if plottype.startswith('lght'):
        posi['ampl'] = np.array([0.1, -0.0])
        posi['sind'] = np.array([0.3, -0.0])
    if plottype == 'lght0002':
        posi['expc'] = np.array([0.5, -0.0])

    if plottype.startswith('lens'):
        posi['defs'] = np.array([0.1, -0.0])
    if plottype == 'lens0001':
        posi['asca'] = np.array([0.3, -0.0])
        posi['acut'] = np.array([0.5, -0.0])
    posi['modl'] = np.array([0., -0.075])
    posi['data'] = np.array([0., -0.15])
   
    if verbtype > 0:
        numb = max(len(grap.edges()), len(listcolr))
        for k in range(numb):
            try:
                print '%15s %15s %15s' % (grap.edges()[k][0], grap.edges()[k][1], listcolr[k])
            except:
                print 'unequal'
        print

    size = 1000
    nx.draw(grap, posi, labels=labl, ax=axis, edgelist=[], nodelist=[])
    nx.draw_networkx_edges(grap, posi, ax=axis, labels=labl, edge_color=listcolr)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['modl', 'data'], node_color='grey', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['numbelem'], node_color='b', node_size=size)
    if plottype.startswith('lght'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'ampldistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'ampl', 'sind'], node_color='g', node_size=size)
    if plottype == 'lght0001' or plottype == 'lght0002':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['sinddistmean'], node_color='r', node_size=size)
    if plottype == 'lght0002':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['expc'], node_color='g', node_size=size)
    if plottype == 'lght0003':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['spatdistcons'], node_color='r', node_size=size)
    nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['psfp', 'bacp'], node_color='y', node_size=size)
    if plottype.startswith('lens'):
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['meanelem', 'defsdistslop'], node_color='r', node_size=size)
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lenp'], node_color='y', node_size=size)
    if plottype == 'lens0000':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs'], node_color='g', node_size=size)
    if plottype == 'lens0001':
        nx.draw_networkx_nodes(grap, posi, ax=axis, labels=labl, nodelist=['lgal', 'bgal', 'defs', 'asca', 'acut'], node_color='g', node_size=size)
    
    pathplot = os.environ["PCAT_DATA_PATH"] + '/imag/'
    plt.tight_layout()
    figr.savefig(pathplot + 'grap%s.pdf' % plottype)
    plt.close(figr)


def plot_3fgl_thrs(gdat):

    path = os.environ["PCAT_DATA_PATH"] + '/detthresh_P7v15source_4years_PL22.fits'
    fluxthrs = astropy.io.fits.getdata(path, 0)

    bgalfgl3 = np.linspace(-90., 90., 481)
    lgalfgl3 = np.linspace(-180., 180., 960)

    bgalexpo = np.linspace(-90., 90., 400)
    lgalexpo = np.linspace(-180., 180., 800)

    #fluxthrs = interp2d(lgalfgl3, bgalfgl3, fluxthrs)(lgalexpo, bgalexpo)
    fluxthrs = griddata([lgalfgl3, bgalfgl3], fluxthrs, [gdat.lgalheal])

    cntsthrs = fluxthrs * gdat.expo

    jbgal = np.where(abs(bgalexpo) < 10.)[0]
    jlgal = np.where(abs(lgalexpo) < 10.)[0]
    extent = [-10, 10, -10, 10]
    
    figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
    axis.set_xlabel(gdat.labllgaltotl)
    axis.set_ylabel(gdat.lablbgaltotl)

    imag = plt.imshow(fluxthrs[np.amin(jbgal):np.amax(jbgal)+1, np.amin(jlghprofi):np.amax(jlghprofi)+1], origin='lower', cmap='Reds', extent=gdat.exttrofi)
    plt.colorbar(imag, fraction=0.05)
    plt.tight_layout()
    figr.savefig(gdat.pathplotrtag + 'thrs.pdf')
    plt.close(figr)
    

def plot_init(gdat):
        
    # make initial plots
    if gdat.makeplot:
        
        if gdat.fittnumbtrap > 0:
            for l in gdat.fittindxpopl:
                if (gdat.fittelemspatevaltype[l] != 'full' and gdat.fittmaxmnumbelempopl[l] > 0) and gdat.numbpixl > 1:
                    if gdat.fittboolelemsbrtdfnc[l]:
                        plot_eval(gdat, l)
                    plot_indxprox(gdat)
        
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if gdat.datatype == 'mock' and gdat.truelensmodltype != 'none':
                figr, axis, path = init_figr(gdat, None, 'post', 'cntpmodlraww', 'this', 'true', i, m, -1)
                imag = retr_imag(gdat, axis, gdat.truecntpmodlraww, 'this', 'true', 'cntpdata', i, m, tdim=True)
                make_cbar(gdat, axis, imag, 0, tick=gdat.tickcntpdata, labl=gdat.lablcntpdata)
                plt.tight_layout()
                figr.savefig(path)
                plt.close(figr)

    if gdat.correxpo:
        gdat.lablnumbpixl = r'$N_{\rm{pix}}$'
        gdat.limtexpo = [gdat.minmexpo, gdat.maxmexpo]
        if gdat.enerbins:
            path = gdat.pathinit + 'expototlmean.pdf'
            tdpy.util.plot_gene(path, gdat.meanener, gdat.expototlmean, scalxdat='logt', scalydat='logt', lablxdat=gdat.lablenertotl, \
                                                                                            lablydat=gdat.lablexpototl, limtydat=gdat.limtexpo)
        
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
                axis.hist(gdat.expototl[i, :], gdat.binsexpo)
                axis.set_xlabel(gdat.lablexpototl)
                axis.set_ylabel(gdat.lablnumbpixl)
                axis.set_xscale('log')
                axis.set_yscale('log')
                plt.tight_layout()
                path = gdat.pathinit + 'histexpoen%02d.pdf' % i
                figr.savefig(path)
                plt.close(figr)
        else:
            figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
            axis.hist(gdat.expototl[:, :].flatten(), gdat.binsexpo)
            axis.set_xlabel(gdat.lablexpototl)
            axis.set_ylabel(gdat.lablnumbpixl)
            axis.set_xscale('log')
            axis.set_yscale('log')
            plt.tight_layout()
            path = gdat.pathinit + 'histexpo.pdf'
            figr.savefig(path)
            plt.close(figr)
            
        if gdat.numbpixl > 1:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    figr, axis, path = init_figr(gdat, None, 'post', 'expo', '', '', i, m, -1)
                    imag = retr_imag(gdat, axis, gdat.expo, '', '', 'expo', i, m)
                    make_cbar(gdat, axis, imag, i)
                    plt.tight_layout()
                    figr.savefig(path)
                    plt.close(figr)
                

def plot_defl(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, \
                                        strgvarb='defl', strgcomp='', indxdefl=None, indxpoplplot=-1, multfact=1., indxenerplot=None, indxevttplot=None):

    if indxdefl is not None:
        strgvarb += 'sing'
    strgvarb = strgvarb + strgcomp
    
    defl = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
    
    defl *= multfact
   
    if indxenerplot is not None:
        defl = defl[indxenerplot, :, indxevttplot, ...]

    if indxdefl is not None:
        defl = defl[..., indxdefl]
        strgvarb += '%04d' % indxdefl
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))

    figr, axis, path = init_figr(gdat, gdatmodi, strgpdfn, strgvarb, strgstat, strgmodl, indxenerplot, indxevttplot, indxpoplplot)
    make_legdmaps(gdat, strgstat, strgmodl, axis)
    draw_frambndr(gdat, axis)
  
    defllgal = defl[:, :, 0]
    deflbgal = defl[:, :, 1]
    fact = 4
    ptch = axis.quiver(gdat.anglfact * gdat.lgalgridcart[::fact, ::fact], gdat.anglfact * gdat.bgalgridcart[::fact, ::fact], \
                       gdat.anglfact * defllgal[::fact, ::fact], gdat.anglfact * deflbgal[::fact, ::fact], scale_units='xy', angles='xy', scale=1)
    supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis)
    plt.subplots_adjust(left=0.2, bottom=0.15, top=0.75, right=0.85)
    plt.subplots_adjust()
    savefigr(gdat, gdatmodi, figr, path)
    plt.close(figr)
    

def plot_genemaps(gdat, gdatmodi, strgstat, strgmodl, strgpdfn, strgvarb, indxenerplot=None, indxevttplot=-1, strgcbar=None, \
                                                                tdim=False, indxpoplplot=-1, strgmome='pmea', intreval=False):
    
    if strgcbar is None:
        strgcbar = strgvarb
  
    # construct the string for the map
    if strgvarb == 'cntpdata':
        strgplot = strgvarb
    else:
        if strgstat == 'post':
            strgtemp = strgmome + strgpdfn
        else:
            strgtemp = ''
        strgplot = strgtemp + strgvarb
    
    figr, axis, path = init_figr(gdat, gdatmodi, strgpdfn, strgplot, strgstat, strgmodl, indxenerplot, indxevttplot, indxpoplplot, intreval=intreval)
   
    maps = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
    
    imag = retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot, indxevttplot, tdim=tdim)
    tick = getattr(gdat, 'tick' + strgcbar) 
    labl = getattr(gdat, 'labl' + strgcbar) 

    make_cbar(gdat, axis, imag, tick=tick, labl=labl)
    make_legdmaps(gdat, strgstat, strgmodl, axis)
    if gdat.suprelem:
        supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxpoplplot)

    # Add two sliders for tweaking the parameters
    if intreval:
        print 'Interactive session began...'
        
        freq_slider = [[] for k in gdat.fittindxpara]
        for k, namepara in enumerate(gdat.fittnamepara):
            if k in gdat.fittindxfixpnumbelem or k in gdat.fittindxfixpmeanelem or k in gdat.fittindxfixpdist:
                continue
            initvalu = gdatmodi.thissamp[k]
            if k in gdat.fittindxfixp:
                labl = getattr(gdat, 'labl' + namepara)
                factplot = gdat.fittfactfixpplot[k]
                minm = getattr(gdat, 'minm' + namepara) * factplot
                maxm = getattr(gdat, 'maxm' + namepara) * factplot
            else:
                factplot = getattr(gdat, 'fact' + namepara[:-12] + 'plot')
                labl = '$%s$' % getattr(gdat, 'labl' + namepara[:-12])
                minm = getattr(gdat, 'minm' + namepara[:-12]) * factplot
                maxm = getattr(gdat, 'maxm' + namepara[:-12]) * factplot
            initvalu *= factplot
            freq_slider_ax = figr.add_axes([0.08, 0.02 * k, 0.1, 0.02])
            freq_slider[k] = Slider(freq_slider_ax, labl, minm, maxm, valinit=initvalu)
            freq_slider[k].label.set_size(10)
            
        def sliders_on_changed(val):
            print 'Slider changed.'
            for k, namepara in enumerate(gdat.fittnamepara):
                if k in gdat.fittindxfixpnumbelem or k in gdat.fittindxfixpmeanelem or k in gdat.fittindxfixpdist:
                    continue
                print namepara
                print freq_slider[k].val
                if k in gdat.fittindxfixp:
                    factplot = gdat.fittfactfixpplot[k]
                else:
                    factplot = getattr(gdat, 'fact' + namepara[:-12] + 'plot')
                gdatmodi.thissamp[k] = freq_slider[k].val / factplot
            print
            proc_samp(gdat, gdatmodi, 'this', 'fitt')
            maps = retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, strgpdfn)
            retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxenerplot, indxevttplot, tdim=tdim, imag=imag)
            for ptch in axis.get_children():
                if isinstance(ptch, mpl.patches.Circle) or isinstance(ptch, mpl.collections.PathCollection):#isinstance(ptch, mpl.lines.Line2D):
                    ptch.remove()
            supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxpoplplot)
            plt.show(block=False)

        for k, namepara in enumerate(gdat.fittnamepara):
            if k in gdat.fittindxfixpnumbelem or k in gdat.fittindxfixpmeanelem or k in gdat.fittindxfixpdist:
                continue
            freq_slider[k].on_changed(sliders_on_changed)

        plt.show()
       
        inpt = raw_input("Press enter to continue...")
        plt.close()
        raise Exception('Interactive session ended...')
    else:
        plt.tight_layout()
        savefigr(gdat, gdatmodi, figr, path)
        plt.close(figr)

