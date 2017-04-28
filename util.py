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
                if intpwdth >= amin(psfntemp[indxanglgood]) and intpwdth <= amax(psfntemp[indxanglgood]):
                    wdthtemp = interp1d_pick(psfntemp[indxanglgood], gdat.binsangl[indxanglgood])(intpwdth)
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


def retr_indxsampcomp(gdat, indxpntsfull, strgmodl):
    
    numbtrapcuml = getattr(gdat, strgmodl + 'numbtrapcuml')
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    indxcomp = getattr(gdat, strgmodl + 'indxcomp')
    spectype = getattr(gdat, strgmodl + 'spectype')
    indxsampcompinit = getattr(gdat, strgmodl + 'indxsampcompinit')
    liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
    
    indxsampcomp = dict()
    for strgcomp in liststrgcomptotl:
        indxsampcomp[strgcomp] = [[] for l in indxpopl]

    indxsampcomp['comp'] = []
    for l in indxpopl:
        indxsamptemp = indxsampcompinit + numbtrapcuml[l] + array(indxpntsfull[l], dtype=int) * numbcomp[l]
        indxsampcomp['lgal'][l] = indxsamptemp
        indxsampcomp['bgal'][l] = indxsamptemp + 1
        if gdat.elemtype == 'lght':
            indxsampcomp['flux'][l] = indxsamptemp + 2
            if gdat.numbener > 1:
                indxsampcomp['sind'][l] = indxsamptemp + 3
                if spectype[l] == 'curv':
                    indxsampcomp['curv'][l] = indxsamptemp + 4
                if spectype[l] == 'expo':
                    indxsampcomp['expo'][l] = indxsamptemp + 4
        if gdat.elemtype == 'lens':
            indxsampcomp['defs'][l] = indxsamptemp + 2
            if gdat.variasca:
                indxsampcomp['asca'][l] = indxsamptemp + 3
            if gdat.variacut:
                indxsampcomp['acut'][l] = indxsamptemp + 4
        if gdat.elemtype == 'clus':
            indxsampcomp['nobj'][l] = indxsamptemp + 2
        indxsampcomp['comp'].append(repeat(indxsamptemp, numbcomp[l]) + tile(indxcomp[l], len(indxpntsfull[l])))
             
    return indxsampcomp


def retr_plotpath(gdat, gdatmodi, strg, strgplot, nameinte=''):
    
    if gdatmodi == None:
        if strg == 'true':
            path = gdat.pathinit + nameinte + strgplot + '.pdf'
        else:
            path = gdat.pathplot + gdat.namesampdist + '/finl/' + nameinte + strgplot + '.pdf'
    else:
        path = gdat.pathplot + gdat.namesampdist + '/fram/' + nameinte + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


def retr_pntsfluxtotl(gdat, lgal, bgal, spec, psfnintp, oaxitype, evalcirc):

    pntsflux = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    indxpixleval = retr_indxpixleval(gdat, lgal, bgal, spec, evalcirc)
    
    numbpnts = lgal.size
    for k in range(numbpnts):
        pntsflux[:, indxpixleval, :] += retr_pntsflux(gdat, lgal, bgal, spec, psfnintp, oaxitype, indxpixleval)
    

def retr_indxpixleval(gdat, lgal, bgal, spec, evalcirc):
    
    if not isscalar(lgal):
        lgal = lgal[0]
        bgal = bgal[0]
    if spec.ndim == 2:
        spec = spec[:, 0]
     
    if evalcirc == 'locl':  
        indxpixlpnts = retr_indxpixl(gdat, bgal, lgal)
        indxfluxproxtemp = digitize(spec[gdat.indxenerfluxdist[0]], gdat.binsprox) - 1

        indxpixleval = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
        if isinstance(indxpixleval, int):
            indxpixleval = gdat.indxpixl
    else:
        indxpixleval = gdat.indxpixl
    
    return indxpixleval


def retr_pntsflux(gdat, lgal, bgal, spec, psfnintp, oaxitype, indxpixleval):
    
    # calculate the distance to all pixels from each point source
    dist = retr_angldistunit(gdat, lgal, bgal, indxpixleval)
   
    # interpolate the PSF onto the pixels
    if oaxitype:
        indxoaxitemp = retr_indxoaxipnts(gdat, lgal, bgal)
        psfntemp = psfnintp[indxoaxitemp](dist)
    else:
        psfntemp = psfnintp(dist)
    
    # scale by the PS spectrum
    pntsflux = spec[:, None, None] * psfntemp
                
    return pntsflux


def retr_mapslght(gdat, bacp, pntsflux, tempindx):
    
    modlflux = pntsflux[tempindx]
    for c in gdat.fittindxback:
        if gdat.fittspecback[c] != None:
            norm = gdat.fittspecback[c] * bacp[gdat.fittindxbacpback[c]]
        else:
            norm = bacp[gdat.fittindxbacpback[c]]
        modlflux += norm[:, None, None] * gdat.fittbackflux[c][tempindx]        

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


def cdfn_igam(xdat, slop, cutf):
    
    cdfn = sp.stats.invgamma.cdf(xdat, slop - 1., scale=cutf)
    
    return cdfn


def icdf_powr(unit, minm, maxm, slop):
    
    para = (unit * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))

    return para


def icdf_igam(xdat, slop, cutf):
    
    icdf = sp.stats.invgamma.ppf(xdat, slop - 1., scale=cutf)
    
    return icdf


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


def pdfn_igam(xdat, slop, cutf):
    
    pdfn = sp.stats.invgamma.pdf(xdat, slop - 1., scale=cutf)
    
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


def cdfn_atan(para, minmpara, maxmpara):
    
    paraunit = (arctan(para) - arctan(minmpara)) / (arctan(maxmpara) - arctan(minmpara))
    
    return paraunit


def icdf_atan(paraunit, minmpara, maxmpara):

    para = tan((arctan(maxmpara) - arctan(minmpara)) * paraunit + arctan(minmpara))
    
    return para


def pdfn_atan(para, minmpara, maxmpara):

    pdfn = 1. / (para**2 + 1.) / (arctan(maxmpara) - arctan(minmpara))
    
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

    return fixpunit


def icdf_fixp(gdat, strgmodl, fixpunit, thisindxfixp):

    scalfixp = getattr(gdat, strgmodl + 'scalfixp')[thisindxfixp]
    
    if scalfixp == 'self' or scalfixp == 'logt' or scalfixp == 'atan':
        
        minmfixp = getattr(gdat, strgmodl + 'minmfixp')[thisindxfixp]
        factfixp = getattr(gdat, strgmodl + 'factfixp')[thisindxfixp]

        if scalfixp == 'self':
            fixp = icdf_self(fixpunit, minmfixp, factfixp)
        elif scalfixp == 'logt':
            fixp = icdf_logt(fixpunit, minmfixp, factfixp)
        elif scalfixp == 'atan':
            maxmfixp = getattr(gdat, strgmodl + 'maxmfixp')[thisindxfixp]
            fixp = icdf_atan(fixpunit, minmfixp, maxmfixp)
    
    elif scalfixp == 'gaus' or scalfixp == 'eerr':
        
        meanfixp = getattr(gdat, strgmodl + 'meanfixp')[thisindxfixp]
        stdvfixp = getattr(gdat, strgmodl + 'stdvfixp')[thisindxfixp]
        
        if scalfixp == 'eerr':
            
            cdfnminmfixp = getattr(gdat, strgmodl + 'cdfnminmfixp')[thisindxfixp]
            cdfndifffixp = getattr(gdat, strgmodl + 'cdfndifffixp')[thisindxfixp]
        
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
        gdatmodi.indxpoplmodi = choice(gdat.fittindxpopl)
    else:
        gdatmodi.indxpoplmodi = thisindxpopl

    if (rand() < gdat.probtran or brth or deth) and (gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] != gdat.fittminmnumbpnts[gdatmodi.indxpoplmodi] or \
                                                     gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] != gdat.fittmaxmnumbpnts[gdatmodi.indxpoplmodi]):
        
        if rand() < gdat.probbrde or brth or deth:
            
            ## births and deaths
            if gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.fittmaxmnumbpnts[gdatmodi.indxpoplmodi] or deth:
                gdatmodi.propdeth = True
            elif gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.fittminmnumbpnts[gdatmodi.indxpoplmodi] or brth:
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
            if gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.fittmaxmnumbpnts[gdatmodi.indxpoplmodi]:
                gdatmodi.propmerg = True
            elif gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] == gdat.fittminmnumbpnts[gdatmodi.indxpoplmodi]:
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
            if gdat.propcomp:
                gdatmodi.thisindxsampfull = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp['comp'])))
            else:
                gdatmodi.thisindxsampfull = gdat.indxfixpprop
            gdatmodi.indxsampmodi = choice(gdatmodi.thisindxsampfull)
            gdatmodi.thisindxproptype = gdat.indxstdppara[gdatmodi.indxsampmodi]
        else:
            gdatmodi.thisindxproptype = gdat.indxproptypewith
        gdatmodi.propwith = True
    
    gdatmodi.proptran = gdatmodi.propbrth or gdatmodi.propdeth or gdatmodi.propsplt or gdatmodi.propmerg

    if gdat.verbtype > 1:
        print 
        print 'retr_thisindxprop()'
        print 'propwith'
        print gdatmodi.propwith
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
        if gdatmodi.propwith and gdat.propwithsing:
            print 'indxsampmodi'
            print gdatmodi.indxsampmodi
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


def retr_probpois(data, modl):
    
    lprb = data * log(modl) - modl - sp.special.gammaln(data + 1)
    
    if False:
        print 'retr_probpois'
        print 'data'
        print data
        print 'modl'
        print modl
        print 'lprb'
        print lprb
        print

    return lprb
    
        
def retr_sampvarb(gdat, indxsampcomp, samp, strgmodl):
    
    indxfixpnumbpnts = getattr(gdat, strgmodl + 'indxfixpnumbpnts')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    
    sampvarb = zeros_like(samp)
    sampvarb[indxfixpnumbpnts] = samp[indxfixpnumbpnts]
    
    for k in indxfixp:
        sampvarb[k] = icdf_fixp(gdat, strgmodl, samp[k], k)
    
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    retr_sampvarbcomp(gdat, strgmodl, indxsampcomp, indxpopl, liststrgcomp, listscalcomp, samp, sampvarb)

    return sampvarb
    

def retr_sampvarbcomp(gdat, strgmodl, indxsampcomp, indxpopl, liststrgcomp, listscalcomp, samp, sampvarb, indxelemfull=slice(None)):
    
    for l in indxpopl:
        for k, strgcomp in enumerate(liststrgcomp[l]):
            
            if False:
                print 'retr_sampvarbcomp'
                print 'listscalcomp'
                print listscalcomp
                print 'l'
                print l
                print 'k'
                print k
                print 'strgcomp'
                print strgcomp
                print 'indxelemfull'
                print indxelemfull
                print 'indxsampcomp[strgcomp][l]'
                print indxsampcomp[strgcomp][l]
                print
                print
                print
                print
                print

            if listscalcomp[l][k] == 'self' or listscalcomp[l][k] == 'dexp' or listscalcomp[l][k] == 'expo' or listscalcomp[l][k] == 'powr':
                minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
                if listscalcomp[l][k] == 'powr':
                    maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
                    distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]
                    sampvarb[indxsampcomp[strgcomp][l][indxelemfull]] = icdf_powr(samp[indxsampcomp[strgcomp][l][indxelemfull]], minm, maxm, distslop)
                else:
                    fact = getattr(gdat, strgmodl + 'fact' + strgcomp)
                    sampvarb[indxsampcomp[strgcomp][l][indxelemfull]] = icdf_self(samp[indxsampcomp[strgcomp][l][indxelemfull]], minm, fact)
            if listscalcomp[l][k] == 'igam':
                distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]
                cutf = getattr(gdat, 'cutf' + strgcomp)
                sampvarb[indxsampcomp[strgcomp][l][indxelemfull]] = icdf_igam(samp[indxsampcomp[strgcomp][l][indxelemfull]], distslop, cutf)
            if listscalcomp[l][k] == 'gaus':
                distmean = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distmean')[l]]
                diststdv = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'diststdv')[l]]
                sampvarb[indxsampcomp[strgcomp][l][indxelemfull]] = icdf_gaus(samp[indxsampcomp[strgcomp][l][indxelemfull]], distmean, diststdv)


def retr_mrkrsize(gdat, sign):
    
    minm = getattr(gdat, 'minm' + gdat.namefeatsign) 
    maxm = getattr(gdat, 'maxm' + gdat.namefeatsign) 
    #mrkrsize = (log(sign) - log(minm)) / (log(maxm) - log(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    mrkrsize = (sqrt(sign) - sqrt(minm)) / (sqrt(maxm) - sqrt(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


def retr_hubbpsfn(gdat):

    gdat.exprpsfp = array([0.05]) / gdat.anglfact
    gdat.exproaxitype = False


def retr_sdsspsfn(gdat):
   
    gdat.exprpsfp = array([0.25 / gdat.anglfact, 1.7e6, 1.9, 0.25 / gdat.anglfact, 2.1e6, 2.])
    gdat.exproaxitype = False


def retr_chanpsfn(gdat):

    gdat.exprpsfp = array([0.35 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.e-1, 2.])
    gdat.exproaxitype = True
   

def retr_sdynpsfn(gdat):

    gdat.exprpsfp = array([0.15 / gdat.anglfact])
    gdat.exproaxitype = False
   

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
    gdatmodi.thisindxsampcomp = deepcopy(gdatmodi.nextindxsampcomp)
    
    if gdat.calcllik and gdat.propwithsing:
        if gdat.elemtype == 'lght':
            gdatmodi.thispntsflux = copy(gdatmodi.nextpntsflux)
        if gdat.elemtype == 'lens':
            gdatmodi.thisdeflelem = copy(gdatmodi.nextdeflelem)


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
    gdat.exprspec[0, 0, :] = fluxchansoft * 0.624e9
    gdat.exprspec[0, 1, :] = fluxchanhard * 0.624e9 / 16.
    # temp
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :]

    # temp
    gdat.exprspec[where(gdat.exprspec < 0.)] = 0.

    gdat.exprsind = -log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :]) / log(gdat.meanener[1] / gdat.meanener[0])
    gdat.exprsind[where(logical_not(isfinite(gdat.exprsind)))[0]] = 2.
    
    # temp
    gdat.exprlgal = tile(gdat.exprlgal, (3, 1)) 
    gdat.exprbgal = tile(gdat.exprbgal, (3, 1)) 
    gdat.exprsind = tile(gdat.exprsind, (3, 1)) 


def retr_fermdata(gdat):
    
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'
    fgl3 = pf.getdata(path)
   
    gdat.exprlgal = deg2rad(fgl3['glon'])
    gdat.exprlgal = ((gdat.exprlgal - pi) % (2. * pi)) - pi
    gdat.exprbgal = deg2rad(fgl3['glat'])
    
    gdat.truenumbpntsfull = gdat.exprlgal.size

    gdat.exprlgal = tile(gdat.exprlgal, (3, 1))
    gdat.exprbgal = tile(gdat.exprbgal, (3, 1))
    
    gdat.exprspec = empty((3, gdat.numbener, gdat.truenumbpntsfull))
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
    gdat.exprsind = tile(gdat.exprsind, (3, 1)) 
    gdat.exprcurv = fgl3['beta']
    gdat.exprexpo = fgl3['Cutoff'] * 1e-3
   
    indxtimevari = where((fgl3timevari < 100.) & (gdat.exprspec[0, gdat.indxenerfluxdist[0], :] > gdat.minmflux))[0]
    
    #indxtimevari = where(gdat.exprspec[0, gdat.indxenerfluxdist[0], :] > gdat.minmflux)[0]
    #gdat.exprlgal = gdat.exprlgal[:, indxtimevari]
    #gdat.exprbgal = gdat.exprbgal[:, indxtimevari]
    #gdat.exprsind = gdat.exprsind[:, indxtimevari]
    #gdat.exprcurv = gdat.exprcurv[:, indxtimevari]
    #gdat.exprexpo = gdat.exprexpo[:, indxtimevari]
    #gdat.exprspec = gdat.exprspec[:, :, indxtimevari]


def retr_rtag(gdat):
    
    rtag = '%d' % (gdat.numbswep)
        
    return rtag


def retr_gaus(gdat, gdatmodi, indxparamodi, stdvpara):
    
    numbparamodi = indxparamodi.size
    stdvtemp = normal(size=numbparamodi) * stdvpara
    if rand() < gdat.fracproprand:
        stdvtemp[choice(indxparamodi)] += randn()
    
    gdatmodi.nextsamp[indxparamodi] = gdatmodi.thissamp[indxparamodi] + stdvtemp
    gdatmodi.nextsamp[indxparamodi] = gdatmodi.nextsamp[indxparamodi] % 1.

       
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
    for k in gdat.fittindxpara:
        if k == gdat.fittnumbfixp:
            print
        if k < gdat.fittnumbfixp:
            name = gdat.fittnamefixp[k]
        else:
            name = ''

        gdatmodi.indxparaprop = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp['comp'])))
        if k in concatenate(gdatmodi.thisindxsampcomp['lgal']):
            print
        try:
            strgboolmodi = '%s' % (k in gdatmodi.indxsampmodi)
        except:
            strgboolmodi = ''
        print '%22s %14.4g %14.4g %14.4g %14.4g %14.4g %14s %10d' % (name, getattr(gdatmodi, 'thissamp')[k], getattr(gdatmodi, 'nextsamp')[k], gdatmodi.thissampvarb[k], \
                                               gdatmodi.nextsampvarb[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k], strgboolmodi, gdat.indxstdppara[k])
    

def show_samplong(gdat, gdatmodi):
    
    for strgmodl in ['fitt', 'true']:
        namepara = getattr(gdat, strgmodl + 'namepara')
        scalpara = getattr(gdat, strgmodl + 'scalpara')
        indxpara = getattr(gdat, strgmodl + 'indxpara')
        if strgmodl == 'fitt':
            indxsampcomp = gdatmodi.thisindxsampcomp
            sampvarb = gdatmodi.thissampvarb
        else:
            indxsampcomp = gdat.trueindxsampcomp
            sampvarb = gdat.truesampvarb

        print 'modl: ' + strgmodl
        print '%20s %20s %15s' % ('namepara', 'sampvarb', 'scalpara')
        for k in indxpara:
            if k in concatenate(indxsampcomp['lgal']):
                print
            print '%20s %20f %15s' % (namepara[k], sampvarb[k], scalpara[k])
    

def rscl_elem(gdat, gdatmodi, indxsampmodi=None):
    
    if indxsampmodi != None:
        if indxsampmodi in gdat.fittindxfixpdist:
            indxpoplmodi = int(gdat.fittnamepara[indxsampmodi][-1])
            strgcompmodi = gdat.fittnamepara[indxsampmodi][:4]
            # temp
            if strgcompmodi == 'gang':
                return
            indxcomp = gdat.fittliststrgcomp[indxpoplmodi].index(strgcompmodi)
            scalcompmodi = gdat.fittlistscalcomp[indxpoplmodi][indxcomp]
            listindxpopltemp = [indxpoplmodi]
        else:
            return
    else:
        listindxpopltemp = gdat.fittindxpopl
    
    for l in listindxpopltemp:
        if indxsampmodi != None:
            liststrgcomptemp = [strgcompmodi]
            listscalcomptemp = [scalcompmodi]
        else:
            liststrgcomptemp = gdat.fittliststrgcomp[l]
            listscalcomptemp = gdat.fittlistscalcomp[l]
        for k, strgcomp in enumerate(liststrgcomptemp):
            comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][l]]
            if listscalcomptemp[k] == 'powr' or listscalcomptemp[k] == 'gaus':
                if listscalcomptemp[k] == 'powr' or listscalcomptemp[k] == 'igam':
                    slop = gdatmodi.nextsampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[l]]
                    if listscalcomptemp[k] == 'powr':
                        minm = getattr(gdat, 'fittminm' + strgcomp)
                        maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                        unit = cdfn_powr(comp, minm, maxm, slop)
                    if listscalcomptemp[k] == 'igam':
                        cutf = getattr(gdat, 'cutf' + strgcomp)
                        unit = cdfn_igam(comp, slop, cutf)
                if listscalcomptemp[k] == 'gaus':
                    distmean = gdatmodi.nextsampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[l]]
                    diststdv = gdatmodi.nextsampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[l]]
                    unit = cdfn_gaus(comp, distmean, diststdv)
                gdatmodi.nextsamp[gdatmodi.thisindxsampcomp[strgcomp][l]] = unit


def retr_prop(gdat, gdatmodi, thisindxpnts=None):
 
    if gdat.verbtype > 1:
        print 'retr_prop()'
    
    gdatmodi.nextsamp = copy(gdatmodi.thissamp)
    gdatmodi.nextsampvarb = copy(gdatmodi.thissampvarb)
    gdatmodi.nextindxpntsfull = deepcopy(gdatmodi.thisindxpntsfull)
  
    if gdat.optiprop:
        while True:
            gdatmodi.nextstdvstdp = copy(gdatmodi.thisstdvstdp)
            gdatmodi.nextstdvstdp[gdatmodi.cntrstdpmodi] += randn() * 1e-4
            if gdatmodi.nextstdvstdp[gdatmodi.cntrstdpmodi] > 0.:
                break
        stdvstdp = gdatmodi.nextstdvstdp
    else:
        stdvstdp = gdat.stdvstdp * gdatmodi.thistmprfactstdv

    if gdatmodi.propwith:
        
        if gdat.propwithsing:
            
            if gdatmodi.indxsampmodi in gdat.fittindxfixp:
                gdatmodi.propfixp = True
            else:
                gdatmodi.propfixp = False
                gdatmodi.indxtrapmodi = gdatmodi.indxsampmodi - gdat.fittindxsampcompinit
                gdatmodi.indxpoplmodi = amin(where(gdatmodi.indxtrapmodi // gdat.fittnumbtrapcumr == 0)[0])
                gdatmodi.numbparapoplinit = gdatmodi.indxtrapmodi - gdat.fittnumbtrapcuml[gdatmodi.indxpoplmodi]
                gdatmodi.indxelemmodi = gdatmodi.numbparapoplinit // gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
                gdatmodi.indxelemfullmodi = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].index(gdatmodi.indxelemmodi)
                gdatmodi.indxcompmodi = gdatmodi.numbparapoplinit % gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
                gdatmodi.thisstrgcomp = gdat.fittliststrgcomp[gdatmodi.indxpoplmodi][gdatmodi.indxcompmodi]
                gdatmodi.thisliststrgcomp = [[] for l in gdat.fittindxpopl]
                gdatmodi.thisliststrgcomp[gdatmodi.indxpoplmodi].append(gdatmodi.thisstrgcomp)
                gdatmodi.thislistscalcomp = [[] for l in gdat.fittindxpopl]
                gdatmodi.thislistscalcomp[gdatmodi.indxpoplmodi].append(gdat.fittlistscalcomp[gdatmodi.indxpoplmodi][gdatmodi.indxcompmodi])
               
                if False:
                    print 'gdatmodi.thissampvarb[gdat.fittindxsampnumbpnts]'
                    print gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts]
                    print 'gdatmodi.nextsampvarb[gdat.fittindxsampnumbpnts]'
                    print gdatmodi.nextsampvarb[gdat.fittindxfixpnumbpnts]
                    print 'gdatmodi.thisindxpntsfull'
                    print gdatmodi.thisindxpntsfull
                    print 'gdatmodi.nextindxpntsfull'
                    print gdatmodi.nextindxpntsfull
                    print 'gdatmodi.indxtrapmodi'
                    print gdatmodi.indxtrapmodi
                    print 'gdatmodi.indxpoplmodi'
                    print gdatmodi.indxpoplmodi
                    print 'gdatmodi.numbparapoplinit'
                    print gdatmodi.numbparapoplinit
                    print 'gdatmodi.indxelemfullmodi'
                    print gdatmodi.indxelemfullmodi
                    print 'gdatmodi.indxcompmodi'
                    print gdatmodi.indxcompmodi
                    print 'gdatmodi.thisliststrgcomp'
                    print gdatmodi.thisliststrgcomp

            ## propose a setp in a parameter
            if gdatmodi.propfixp:
                stdvpara = stdvstdp[gdat.indxstdppara[gdatmodi.indxsampmodi]]
            else:
                stdvpara = retr_propcompscal(gdat, gdatmodi, stdvstdp, gdatmodi.indxpoplmodi, gdatmodi.thisstrgcomp, indxelemfull=gdatmodi.indxelemfullmodi)
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, stdvpara)
            
            if gdatmodi.propfixp:
                gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_fixp(gdat, 'fitt', gdatmodi.nextsamp[gdatmodi.indxsampmodi], gdatmodi.indxsampmodi)
                
            if gdatmodi.indxsampmodi in gdat.fittindxfixpdist:
                gdatmodi.propdist = True
                ### rescale the element components due to the hyperparameter step
                rscl_elem(gdat, gdatmodi, gdatmodi.indxsampmodi)
            else:
                gdatmodi.propdist = False
            
            if not gdatmodi.propfixp:
                indxpopl = [gdatmodi.indxpoplmodi]
                
                retr_sampvarbcomp(gdat, 'fitt', gdatmodi.thisindxsampcomp, indxpopl, gdatmodi.thisliststrgcomp, gdatmodi.thislistscalcomp, gdatmodi.nextsamp, \
                                                                                                                gdatmodi.nextsampvarb, indxelemfull=gdatmodi.indxelemfullmodi)
                
            ### asymmetric proposal acceptance ratio factor
            
        else:
            
            ## propose steps in all fixed dimensional, floating parameters
            retr_gaus(gdat, gdatmodi, gdat.indxfixpprop, stdvstdp[gdat.indxstdppara[gdat.indxfixpprop]])
            
            ### rescale the element components due to the hyperparameter steps
            if gdat.propdist:
                for k in gdat.fittindxfixpdist:
                    gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, 'fitt', gdatmodi.nextsamp[k], k)
                rscl_elem(gdat, gdatmodi)

            ### element
            gdatmodi.thislfctprop = 0.
            if gdat.propcomp:
                for l in gdat.fittindxpopl:
                    for strgcomp in gdat.fittliststrgcomp[l]:
                        stdvpara = retr_propcompscal(gdat, gdatmodi, stdvstdp, l, strgcomp)
                        retr_gaus(gdat, gdatmodi, gdatmodi.thisindxsampcomp[strgcomp][l], stdvpara)

    else:
        # number of unit sample vector elements to be modified
        numbcompmodi = gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
    
    if gdatmodi.propbrth or gdatmodi.propsplt:
       
        if gdatmodi.propbrth:
            numbiter = 1
        else:
            numbiter = 2
    
        gdatmodi.indxsamptran = []
        cntr = 0
        for a in range(numbiter):
            # find an empty slot in the PS list
            for k in range(gdat.fittmaxmnumbpnts[gdatmodi.indxpoplmodi]):
                if not k in gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi]:
                    indxpntsbrth = k
                    cntr += 1
                if cntr == numbiter:
                    break
       
            # sample indices to add the new PS
            gdatmodi.indxsamptran.append(retr_indxsamppnts(gdat, gdatmodi.indxpoplmodi, array([indxpntsbrth])))
        
        # sample auxiliary variables
        gdatmodi.auxipara = rand(numbcompmodi)
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.auxipara
                
        # change the number of PS
        gdatmodi.nextsamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].append(indxpntsbrth)
    
    # death
    if gdatmodi.propdeth:
        
        if thisindxpnts != None:
            dethindxindxpnts = thisindxpnts
        else:
            # occupied PS index to be killed
            dethindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
            
        # PS index to be killed
        indxpntsdeth = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][dethindxindxpnts]
    
        # sample indices to add the new PS
        gdatmodi.indxsamptran = []
        gdatmodi.indxsamptran.append(gdat.fittindxsampcompinit + gdat.fittnumbtrapcuml[gdatmodi.indxpoplmodi] + indxpntsdeth * gdat.fittnumbcomp[gdatmodi.indxpoplmodi] + \
                                                                                                                                    gdat.fittindxcomp[gdatmodi.indxpoplmodi])
    
        # remove the PS from the occupied PS list
        gdatmodi.nextindxpntsfull[gdatmodi.indxpoplmodi].remove(indxpntsdeth)
    
        # change the number of PS
        gdatmodi.nextsamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
                  
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbpntsmodi = 3
        
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        gdatmodi.spltindxindxpnts = choice(thisindxindxpnts)
    
        # update the sample vector
        gdatmodi.indxsampfrst = gdat.fittindxsampcompinit + gdat.fittnumbtrap * gdatmodi.indxpoplmodi + \
                                                int(gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]) * gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        indxfinlfrst = gdatmodi.indxsampfrst + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]

        gdatmodi.indxsampseco = gdat.fittindxsampcompinit + gdat.fittnumbtrap * gdatmodi.indxpoplmodi + indxpntsbrth * gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        indxfinlseco = gdatmodi.indxsampseco + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, indxfinlfrst, dtype=int), arange(gdatmodi.indxsampseco, indxfinlseco, dtype=int)))
        
        for k, strgcomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
            if strgcomp == 'lgal':
                gdatmodi.lgalpare = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
            if strgcomp == 'bgal':
                gdatmodi.bgalpare = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
            if strgcomp == gdat.namefeatsign:
                gdatmodi.bgalpare = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
                setattr(gdatmodi, gdatmodi.bgalpare, compsign)
            
        gdatmodi.fluxpare = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[''][gdatmodi.indxpoplmodi][gdatmodi.spltindxindxpnts]]
        
        # determine the new element parameters
        # temp -- only valid for power-law energy spectrum
        gdatmodi.auxipara = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.auxipara[0] = rand()
        gdatmodi.auxipara[1] = rand() * gdat.radispmr
        gdatmodi.auxipara[2] = rand() * gdat.radispmr
        # temp
        if gdat.numbener > 1:
            gdatmodi.auxipara[3] = icdf_gaus(rand(), gdat.sinddistmean[gdatmodi.indxpoplmodi], gdat.sinddiststdv[gdatmodi.indxpoplmodi])
            
        gdatmodi.fluxfrst = gdatmodi.auxipara[0] * gdatmodi.fluxpare
        gdatmodi.spltlgalfrst = thislgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1]
        gdatmodi.spltbgalfrst = thisbgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[2]
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
        gdatmodi.thislistpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], \
                                                                                    gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]])
        gdatmodi.thisnumbpair = len(gdatmodi.thislistpair)
        
        if gdatmodi.thisaccpprio:

            # calculate the list of pairs
            ## proposed
            lgal = concatenate((array([gdatmodi.spltlgalfrst, gdatmodi.spltlgalseco]), \
                                                                    setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([gdatmodi.spltbgalfrst, gdatmodi.spltbgalseco]), \
                                                                    setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]], thisbgal)))
            gdatmodi.nextlistpair = retr_listpair(gdat, lgal, bgal)
            gdatmodi.nextnumbpair = len(gdatmodi.nextlistpair)

            if gdatmodi.nextnumbpair == 0:
                raise Exception('Number of pairs should not be zero in the reverse proposal of a split')

            ## first new element
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcomplgal] = cdfn_self(gdatmodi.spltlgalfrst, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompbgal] = cdfn_self(gdatmodi.spltbgalfrst, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompflux] = cdfn_flux_powr(gdatmodi.fluxfrst, gdat.minmflux, gdat.maxmflux, \
                                                                                gdatmodi.thissampvarb[gdat.fittindxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.nextsamp[gdatmodi.indxsampfrst+gdat.indxcompsind] = cdfn_gaus(gdatmodi.spltsindfrst, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
    
            nextspecfrst = retr_spec(gdat, gdatmodi.fluxfrst, spep=gdatmodi.spltsindfrst, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

            ## second new element
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcomplgal] = cdfn_self(gdatmodi.spltlgalseco, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcompbgal] = cdfn_self(gdatmodi.spltbgalseco, -gdat.maxmgang, 2. * gdat.maxmgang)
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcompflux] = cdfn_flux_powr(gdatmodi.fluxseco, gdat.minmflux, gdat.maxmflux, \
                                                                                gdatmodi.thissampvarb[gdat.fittindxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.nextsamp[gdatmodi.indxsampseco+gdat.indxcompsind] = cdfn_gaus(gdatmodi.spltsindseco, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                                                                                                            gdat.sinddiststdv[gdatmodi.indxpoplmodi])
            nextspecseco = retr_spec(gdat, gdatmodi.fluxseco, spep=gdatmodi.spltsindseco, spectype=gdat.spectype[gdatmodi.indxpoplmodi])

    if gdatmodi.propmerg:
        
        # number of point sources to be modified
        gdatmodi.numbpntsmodi = 3
        
        # proposed number of point sources
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # calculate the current list of pairs
        gdatmodi.thislistpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], \
                                                                                        gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]])
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
            gdatmodi.indxsampfrst = gdat.fittindxsampcompinit + gdat.fittnumbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxfrst
            indxfinlfrst = gdatmodi.indxsampfrst + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
            
            ## second PS
            gdatmodi.indxsampseco = gdat.fittindxsampcompinit + gdat.fittnumbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxseco
            indxfinlseco = gdatmodi.indxsampseco + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]

            # indices of the sample vector elements to be modified
            gdatmodi.indxsampmodi = arange(gdatmodi.indxsampfrst, indxfinlfrst)

            # indices of the PS to be merged
            mergindxpnts = sort(array([gdatmodi.mergindxfrst, gdatmodi.mergindxseco], dtype=int))

            # PS parameters to be merged
            ## first PS
            gdatmodi.lgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.bgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            gdatmodi.fluxfrst = gdatmodi.specfrst[gdat.indxenerfluxdist[0]]

            ## second PS
            gdatmodi.lgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.bgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            gdatmodi.fluxseco = gdatmodi.specseco[gdat.indxenerfluxdist[0]]

            # auxiliary parameters
            auxifrac = gdatmodi.fluxfrst / (gdatmodi.fluxfrst + gdatmodi.fluxseco) 
            auxiradi = sqrt((gdatmodi.lgalseco - gdatmodi.lgalfrst)**2 + (gdatmodi.bgalseco - gdatmodi.bgalfrst)**2)
            auxiangl = pi + arctan2(gdatmodi.bgalseco - gdatmodi.bgalfrst, gdatmodi.lgalseco - gdatmodi.lgalfrst)
            auxispep = gdatmodi.spepseco

            # temp
            gdatmodi.auxipara = zeros(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
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
                                                                                        gdatmodi.thissampvarb[gdat.fittindxfixpfluxdistslop[gdatmodi.indxpoplmodi]])
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
            
            lgal = concatenate((array([gdatmodi.lgalpare]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], \
                                                                                                        array([gdatmodi.lgalfrst, gdatmodi.lgalseco]))))
            bgal = concatenate((array([gdatmodi.bgalpare]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]], \
                                                                                                        array([gdatmodi.bgalfrst, gdatmodi.bgalseco]))))
            gdatmodi.nextlistpair = retr_listpair(gdat, lgal, bgal)
            gdatmodi.nextnumbpair = len(gdatmodi.nextlistpair)
        
    if gdat.fittnumbtrap > 0:
        if gdatmodi.propwith:
            if gdat.propwithsing:
                if gdatmodi.propdist:
                    gdatmodi.indxsampmodi = array([gdatmodi.indxsampmodi, ])
            else:
                gdatmodi.indxsampmodi = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp['comp'])))
            gdatmodi.indxsampchec = gdatmodi.indxsampmodi
        else:
            gdatmodi.indxsampchec = []
            if gdatmodi.propbrth:
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[0]))
            if gdatmodi.propdeth:
                gdatmodi.indxsampmodi = gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None]
            if gdatmodi.propsplt:
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[1], gdatmodi.indxsamptran[2]))
            if gdatmodi.propmerg:
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[0]))
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
    
    if gdatmodi.propwith:
        gdatmodi.nextindxsampcomp = gdatmodi.thisindxsampcomp
    else:
        gdatmodi.nextindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.nextindxpntsfull, 'fitt')

    if gdatmodi.thisaccpprio and (gdatmodi.propbrth or gdatmodi.propwith and not gdat.propwithsing):
        gdatmodi.nextsampvarb = retr_sampvarb(gdat, gdatmodi.nextindxsampcomp, gdatmodi.nextsamp, 'fitt')
   
    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        
        ## Jacobian
        jcbnfacttemp = log(gdatmodi.fluxpare * fabs(gdatmodi.auxipara[1] * (sin(gdatmodi.auxipara[2]) * cos(gdatmodi.auxipara[2]) + cos(gdatmodi.auxipara[2])**2)))
        if gdatmodi.propsplt:
            gdatmodi.thisjcbnfact = jcbnfacttemp
        else:
            gdatmodi.thisjcbnfact = -jcbnfacttemp
        
        ## combinatorial factor
        thisnumbpnts = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]]
        if gdatmodi.propsplt:
            gdatmodi.thiscombfact = log(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2 / gdatmodi.nextnumbpair)
        else:
            gdatmodi.thiscombfact = log(gdatmodi.thisnumbpair / gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2)
   

def retr_propcompscal(gdat, gdatmodi, stdvstdp, l, strgcomp, indxelemfull=slice(None)):
    
    thiscompsign = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namefeatsign][l][indxelemfull]]
    nextcompsign = gdatmodi.nextsampvarb[gdatmodi.thisindxsampcomp[gdat.namefeatsign][l][indxelemfull]]
    minmcompsign = getattr(gdat, 'minm' + gdat.namefeatsign)
    stdvstdpcomp = stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)]
    thiscompunit = gdatmodi.thissamp[gdatmodi.thisindxsampcomp[gdat.namefeatsign][l][indxelemfull]]
    nextcompunit = gdatmodi.nextsamp[gdatmodi.thisindxsampcomp[gdat.namefeatsign][l][indxelemfull]]
    if strgcomp == gdat.namefeatsign:
        # temp -- this only works if compsign is powr distributed
        gdatmodi.thisstdv = stdvstdpcomp / (thiscompsign / minmcompsign)**2.
        gdatmodi.nextstdv = stdvstdpcomp / (nextcompsign / minmcompsign)**2.
        gdatmodi.thislfctprop += sum(0.5 * (nextcompunit - thiscompunit)**2 * (1. / gdatmodi.thisstdv**2 - 1. / gdatmodi.nextstdv**2))
    else:
        gdatmodi.thisstdv = stdvstdpcomp / (minimum(thiscompsign, nextcompsign) / minmcompsign)**0.5

    return gdatmodi.thisstdv


def show_dbug(gdat, gdatmodi):
    print 'hey'
    print 'thissamp nextsamp diffsampvarb'
    for k in gdat.fittindxpara:
        print '%10.3g  %10.3g %10.3g' % (gdatmodi.thissamp[k], gdatmodi.nextsamp[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k])
    print

        
def retr_indxsamppnts(gdat, l, indxpnts):

    indxsamppnts = gdat.fittindxsampcompinit + gdat.fittnumbtrapcuml[l] + indxpnts[None, :] * gdat.fittnumbcomp[l] + gdat.fittindxcomp[l][:, None]
    indxsamppnts = indxsamppnts.flatten()

    return indxsamppnts


def retr_factoaxi(gdat, bins, norm, indx):

    factoaxi = 1. + norm[:, None, None] * (bins[None, None, :] / gdat.oaxipivt)**indx[:, None, None]
    
    return factoaxi


def gang_detr():

    gang, aang, lgal, bgal = sympy.symbols('gang aang lgal bgal')

    AB = sympy.matrices.Matrix([[a1*b1,a1*b2,a1*b3],[a2*b1,a2*b2,a2*b3],[a3*b1,a3*b2,a3*b3]])
    print AB.det()


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, binsoaxi=None, oaxitype=None, strgmodl='fitt'):

    numbpsfpform = getattr(gdat, strgmodl + 'numbpsfpform')
    numbpsfptotl = getattr(gdat, strgmodl + 'numbpsfptotl')
    
    indxpsfpinit = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    if oaxitype:
        indxpsfponor = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp]
        indxpsfpoind = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp] + 1
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * arcsin(sqrt(2. - 2. * cos(gdat.binsangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
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
        fact = 2. * pi * trapz(psfnnorm * sin(gdat.binsangl[None, :, None]), gdat.binsangl, axis=1)[:, None, :]
        psfn /= fact

    # temp
    if True and (gdat.strgcnfg == 'pcat_ferm_expr_ngal' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp1' or \
                                                            gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp2' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp3'):
        print 'CORRECTING THE PSF.'
        tempcorr = array([1., 0.8, 0.8])
        psfn *= tempcorr[:, None, None]

    return psfn


def retr_axis(gdat, strgvarb, minm=None, maxm=None, numb=None, bins=None, scal='self', invr=False, strginit=''):
    
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

    if scal == 'logt':
        meanvarb = sqrt(bins[1:] * bins[:-1])
    else:
        meanvarb = (bins[1:] + bins[:-1]) / 2.
    
    indx = arange(numb)
    delt = diff(bins) 

    setattr(gdat, strginit + 'bins' + strgvarb, bins)
    setattr(gdat, strginit + 'mean' + strgvarb, meanvarb)
    setattr(gdat, strginit + 'delt' + strgvarb, delt)
    setattr(gdat, strginit + 'numb' + strgvarb, numb)
    setattr(gdat, strginit + 'indx' + strgvarb, indx)


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
    for k in gdat.indxwvecodim:
        indxwvec = where((gdat.meanwvec > gdat.binswvecodim[k]) & (gdat.meanwvec < gdat.binswvecodim[k+1]))
        psecodim[k] = mean(psec[indxwvec])
    
    psecodim *= gdat.meanwvecodim**2

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


def retr_varb(gdat, gdatmodi, strg, strgvarb, indx=None, perc='medi'):
        
    if strg == 'this':
        #varb = gdatmodi.thissampvarb[getattr(gdat, 'indx' + strgvarb)]
        varb = getattr(gdatmodi, strg + strgvarb)
    elif strg == 'true':
        varb = getattr(gdat, strg + strgvarb)
    else:
        varb = getattr(gdat, perc + strgvarb)
        
    if indx != None:
        if perc == 'errr':
            varb = varb[[slice(None)] + indx]
        else:
            varb = varb[indx]
    
    return varb


def retr_eerrnorm(minmvarb, maxmvarb, meanvarb, stdvvarb):
   
    cdfnminm = 0.5 * (sp.special.erf((minmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfnmaxm = 0.5 * (sp.special.erf((maxmvarb - meanvarb) / stdvvarb / sqrt(2.)) + 1.)
    cdfndiff = cdfnmaxm - cdfnminm
    
    return cdfnminm, cdfndiff
    

def retr_condcatl(gdat):
  
    if gdat.verbtype > 1:
        print 'gdat.listindxpntsfull'
        print gdat.listindxpntsfull

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
            for k in range(len(gdat.listindxpntsfull[n][l])):
                indxstks[n][l].append(numbstks)
                indxstkssamptemp.append(numbstks)
                indxtupl.append([n, l, k])
                numbstks += 1
        indxstkssamp.append(array(indxstkssamptemp))
    
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
    arrystks = zeros((numbstks, gdat.fittnumbcomptotl))
    for n in gdat.indxsamptotl:
        indxsampcomp = retr_indxsampcomp(gdat, gdat.listindxpntsfull[n], 'fitt') 
        for l in gdat.fittindxpopl:
            for k in arange(len(gdat.listindxpntsfull[n][l])):
                for m, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                    arrystks[indxstks[n][l][k], m] = gdat.listsamp[n, indxsampcomp[strgcomp][l][k]]

    if gdat.verbtype > 0:
        print 'Constructing the distance matrix for %d stacked samples...' % arrystks.shape[0]
        timeinit = gdat.functime()
    
    gdat.distthrs = empty(gdat.fittnumbcomptotl)
    for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
        gdat.distthrs[k] = gdat.stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)]
    
    # construct lists of samples for each proposal type

    listdisttemp = [[] for k in range(numbstks)]
    indxstksrows = [[] for k in range(numbstks)]
    indxstkscols = [[] for k in range(numbstks)]
    thisperc = 0
    for k in gdat.fittindxcomptotl:
        for n in range(numbstks):
            dist = fabs(arrystks[n, k] - arrystks[:, k])
            indxstks = where(dist < gdat.distthrs[k])[0]
            if indxstks.size > 0:
                for j in indxstks:
                    listdisttemp[k].append(dist[j])
                    indxstksrows[k].append(n)
                    indxstkscols[k].append(j)
            
            nextperc = floor(100. * float(k * numbstks + n) / numbstks / gdat.fittnumbcomptotl)
            if nextperc > thisperc:
                thisperc = nextperc
                print '%3d%% completed.' % thisperc
        
        listdisttemp[k] = array(listdisttemp[k])
        indxstksrows[k] = array(indxstksrows[k])
        indxstkscols[k] = array(indxstkscols[k])
    
    listdist = [[] for k in range(numbstks)]
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
        numbdist = zeros(numbstks, dtype=int) - 1
        for p in range(len(indxstksleft)):
            indxstksdist = arange(numbstks)
            for k in gdat.fittindxcomptotl:
                indxindx = where(listdist[k][indxstksleft[p], :] < gdat.distthrs[k])[0]
                indxstksdist = intersect1d(indxstksdist, indxstkscols[k][indxindx])
            numbdist[indxstksleft[p]] = indxstksdist.size
            
        prvlmaxmesti = amax(numbdist) / float(gdat.numbsamptotl)
        
        if True or cntr % 10 == 0:
            print 'cntr'
            print cntr
            print 'numbdist'
            print numbdist
            print 'gdat.numbsamptotl'
            print gdat.numbsamptotl
            print 'prvlmaxmesti'
            print prvlmaxmesti

        if prvlmaxmesti < gdat.prvlthrs:
            print 'Exiting...'
            break

        # determine the element with the highest number of neighbors
        indxstkscntr = argmax(numbdist)
        indxsamptotlcntr = indxtupl[indxstkscntr][0]
        indxpoplcntr = indxtupl[indxstkscntr][1]
        indxpntscntr = indxtupl[indxstkscntr][2]

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
                
                indxstkstemp = intersect1d(array(indxstksleft), indxstkssamp[n])
                
                if gdat.verbtype > 1:
                    print 'n'
                    print n
                    print 'indxstkstemp'
                    print indxstkstemp
                
                if n == indxsamptotlcntr:
                    continue
                
                if indxstkstemp.size > 0:
                    totl = zeros_like(indxstkstemp)
                    for k in gdat.fittindxcomptotl:
                        temp = listdist[k][indxstkscntr, indxstkstemp].toarray()[0]
                        totl = totl + temp**2

                    indxleft = argsort(totl)[0]
                    
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
                #    gdatmodi.thisindxpntsfull = deepcopy(gdat.listindxpntsfull[n])
                #    for r in range(len(indxstksassc)): 
                #        calc_poststkscond(gdat, indxstksassc)
                #    gdatmodi.thisindxpntsfull = [[] for l in gdat.fittindxpopl]
                #    for indxstkstemp in indxstksleft:
                #        indxsamptotlcntr = indxtupl[indxstkstemp][0]
                #        indxpoplcntr = indxtupl[indxstkstemp][1]
                #        indxpntscntr = indxtupl[indxstkstemp][2]
                #        gdatmodi.thissampvarb = gdat.listsampvarb[indxsamptotlcntr, :]
                #        gdatmodi.thisindxpntsfull[].append()

                #    plot_genemaps(gdat, gdatmodi, 'this', 'datacnts', thisindxener=0, thisindxevtt=0, cond=True)
                
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
    for r in range(len(indxstksassc)): 
        gdat.dictglob['liststkscond'].append([])
        gdat.dictglob['liststkscond'][r] = {}
        gdat.dictglob['poststkscond'].append([])
        gdat.dictglob['poststkscond'][r] = {}
        for strgfeat in gdat.fittliststrgfeatodimtotl:
            gdat.dictglob['liststkscond'][r][strgfeat] = []
            gdat.dictglob['poststkscond'][r][strgfeat] = zeros(3)
        for k in range(len(indxstksassc[r])):
            indxsamptotlcntr = indxtupl[indxstksassc[r][k]][0]
            indxpoplcntr = indxtupl[indxstksassc[r][k]][1]
            indxpntscntr = indxtupl[indxstksassc[r][k]][2]
            
            if gdat.verbtype > 1:
                print 'indxsamptotlcntr'
                print indxsamptotlcntr
                print 'indxpoplcntr'
                print indxpoplcntr
                print 'indxpntscntr'
                print indxpntscntr
                print
            
            for strgfeat in gdat.fittliststrgfeatodimtotl:
                temp = getattr(gdat, 'list' + strgfeat)
                temp = temp[indxsamptotlcntr][indxpoplcntr][indxpntscntr]
                gdat.dictglob['liststkscond'][r][strgfeat].append(temp)
    
    for r in range(len(gdat.dictglob['liststkscond'])): 
        for strgfeat in gdat.fittliststrgfeatodimtotl:
            arry = array(gdat.dictglob['liststkscond'][r][strgfeat])
            gdat.dictglob['poststkscond'][r][strgfeat][0] = median(arry)
            gdat.dictglob['poststkscond'][r][strgfeat][1] = percentile(arry, 16.)
            gdat.dictglob['poststkscond'][r][strgfeat][2] = percentile(arry, 84.)
    

    gdat.numbstkscond = len(gdat.dictglob['liststkscond'])
    print 'gdat.numbstkscond'
    print gdat.numbstkscond
    print

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
    
    gdat.indxstkscond = arange(gdat.numbstkscond)
    gdat.prvl = empty(gdat.numbstkscond)
    for r in gdat.indxstkscond:
        gdat.prvl[r] = len(gdat.dictglob['liststkscond'][r]['deltllik'])
    gdat.prvl /= gdat.numbsamptotl
    print 'gdat.prvl'
    print gdat.prvl
    gdat.minmprvl = 0.
    gdat.maxmprvl = 1.
    retr_axis(gdat, 'prvl', gdat.minmprvl, gdat.maxmprvl, 10)
    gdat.histprvl = histogram(gdat.prvl, bins=gdat.binsprvl)[0]
    if gdat.makeplot:
        pathcond = getattr(gdat, 'path' + gdat.namesampdist + 'finlcond')
        for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
            path = pathcond + 'histdist' + strgcomp 
            listtemp = copy(listdist[k]).flatten()
            listtemp = listtemp[where(listtemp != 1e20)[0]]
            tdpy.mcmc.plot_hist(path, listtemp, r'$Delta$ ' + getattr(gdat, 'labl' + strgcomp))
            path = pathcond + 'histprvl'
            tdpy.mcmc.plot_hist(path, gdat.prvl, r'$p$')
    gdat.prvlthrs = 0.01 
    gdat.indxprvlhigh = where(gdat.prvl > gdat.prvlthrs)[0]
    gdat.numbprvlhigh = gdat.indxprvlhigh.size
    print 'gdat.numbprvlhigh'
    print gdat.numbprvlhigh


def retr_conv(gdat, defl):
    
    # temp
    conv = abs(gradient(defl[:, :, 0], gdat.sizepixl, axis=0) + gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) / 2.

    return conv


def retr_invm(gdat, defl):
    
    # temp
    invm = (1. - gradient(defl[:, :, 0], gdat.sizepixl, axis=0)) * (1. - gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) - \
                                                                         gradient(defl[:, :, 0], gdat.sizepixl, axis=1) * gradient(defl[:, :, 1], gdat.sizepixl, axis=0)

    return invm


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


def retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut):
    
    critmden = retr_critmden(gdat, adissour, adishost, adishostsour)
    
    fracacutasca = acut / asca
    factmcutfromdefs = pi * adishost**2 * critmden * asca / (1. - log(2.)) * retr_mcutfrommscl(fracacutasca)
        
    return factmcutfromdefs


def retr_mcut(gdat, defs, asca, acut):
    
    mscl = defs * pi * gdat.adishost**2 * gdat.critmden * asca / (1. - log(2.))
    fracacutasca = acut / asca
    mcut = mscl * retr_mcutfrommscl(fracacutasca)
    
    return mcut


def retr_mcutfrommscl(fracacutasca):
    
    mcut = fracacutasca**2 / (fracacutasca**2 + 1.)**2 * ((fracacutasca**2 - 1.) * log(fracacutasca) + fracacutasca * pi - (fracacutasca**2 + 1.))

    return mcut


def setpinit(gdat, boolinitsetp=False):

    if False and gdat.elemtype == 'lens' and gdat.strgproc == 'fink2.rc.fas.harvard.edu':
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
    
    gdat.liststrgmodl = ['fitt', 'true']
    
    gdat.liststrgtype = ['spatdisttype', 'fluxdisttype', 'spectype']
        
    # paths
    ## data
    gdat.pathpixlcnvt = gdat.pathdata + 'pixlcnvt/'
    gdat.pathprox = gdat.pathdata + 'prox/'
    ## plot
    gdat.pathplot = gdat.pathimag + gdat.strgtimestmp + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
    gdat.pathinit = gdat.pathplot + 'init/'
    gdat.pathinitintr = gdat.pathinit + 'intr/'
    
    gdat.pathprio = gdat.pathplot + 'prio/'
    gdat.pathpost = gdat.pathplot + 'post/'
    
    for name in ['prio', 'post']:
        
        path = getattr(gdat, 'path' + name)
        for nameseco in ['finl', 'fram', 'anim', 'opti']:
            setattr(gdat, 'path' + name + nameseco, path + nameseco + '/')
        
        for nameseco in ['diag', 'lpri', 'varbscal', 'cond', 'varbscalproc', 'deltllik', 'spmr']:
            setattr(gdat, 'path' + name + 'finl' + nameseco, path + 'finl/' + nameseco + '/')
            
        for nameseco in ['histodim', 'histtdim', 'assc', 'scattdim']:
            for namethrd in ['init', 'fram', 'finl', 'anim']:
                if namethrd == 'init':
                    if nameseco == 'assc' or nameseco == 'histtdim':
                        continue
                    setattr(gdat, 'path' + namethrd + nameseco, gdat.pathplot + 'init/' + nameseco + '/')
                else:
                    setattr(gdat, 'path' + name + namethrd + nameseco, path + namethrd + '/' + nameseco + '/')
            
    gdat.pathinfo = gdat.pathplot + 'info/'
    
    ## make the directories 
    for attr, valu in gdat.__dict__.iteritems():
        if attr.startswith('path'):
            os.system('mkdir -p %s' % valu)
 
    if gdat.elemtype == 'lens':
        gdat.ascaglob = 0.2 / gdat.anglfact
        gdat.acutglob = 0.6 / gdat.anglfact
    
    gdat.cutfdefs = 3e-3 / gdat.anglfact

    # plotting
    gdat.legdsampdist = 'Posterior'
    gdat.legdsamp = 'Sample'
    gdat.legdmedi = 'Median'
    gdat.legdstdv = 'Std. dev.'
    
    # optimization period
    gdat.numbswepoptiprop = 10 * gdat.fittnumbpara

    # p value contours 
    gdat.pvalcont = [0.317, 0.0455, 2.7e-3, 6e-5, 1.3e-6]

    ## number of bins in histogram plots
    gdat.numbbinsplot = 20
    ## number of bins in hyperprior plots
    gdat.numbbinsplotprio = 100
    # temp
    ### number of standard deviations away from the mean of spectral parameters to plot
    gdat.numbstdvspepdist = 3.
    ### minima and maxima for spectral parameters
    gdat.numbstdv = 3.
    
    ## labels and scales for variables
    if gdat.elemtype == 'lens':
        gdat.lablfracsubh = '$f_{sub}$'

        gdat.lablmasssubhtotl = '$M_{sub}$'
        gdat.scalmasssubhtotl = 'logt'
        
        gdat.lablmasshostbein = '$M_{E,hst}$'
        gdat.scalmasshosttotl = 'logt'
    
    gdat.scalmaxmnumbpnts = 'logt'
    gdat.scalmedilliktotl = 'logt'

    gdat.lablener = 'E'
    #gdat.lablenertotl = '$%s$ [%s]' % (gdat.lablener, gdat.strgenerunit)
    
    gdat.lablgang = r'\theta'
    gdat.lablaang = r'\phi'
    gdat.labllgalunit = gdat.lablgangunit
    gdat.lablbgalunit = gdat.lablgangunit
    gdat.lablaangunit = ''
    
    gdat.labldefs = r'\alpha_s'
    gdat.lablflux = 'f'
    gdat.lablnobj = 'p'
    
    gdat.labldeflprof = r'\alpha_a'
    gdat.labldeflprofunit = u'$^{\prime\prime}$'
    
    gdat.labldefsunit = u'$^{\prime\prime}$'
    if gdat.exprtype == 'ferm':
        gdat.lablfluxunit = 'cm$^{-2}$ s$^{-1}$ GeV$^{-1}$'
        gdat.lablfluxsoldunit = 'cm$^{-2}$ s$^{-1}$ GeV$^{-1} sr$^{-1}$$'
    if gdat.exprtype == 'chan':
        gdat.lablfluxunit = 'cm$^{-2}$ s$^{-1}$ KeV$^{-1}$'
        gdat.lablfluxsoldunit = 'cm$^{-2}$ s$^{-1}$ KeV$^{-1} sr$^{-1}$$'
    if gdat.exprtype == 'hubb':
        gdat.lablfluxunit = r'erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$'
        gdat.lablfluxsoldunit = r'erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$ sr$^{-1}$'
    
    gdat.lablfluxsold = 'f'
    gdat.lablnobjunit = ''
    
    for l in gdat.trueindxpopl:
        setattr(gdat, 'lablcmplpop%d' % l, '$c_{%d}$' % l)
        setattr(gdat, 'lablfdispop%d' % l, '$f_{%d}$' % l)

    gdat.lablprvl = '$p$'
    
    gdat.lablsind = 's'
    gdat.lablsindunit = ''
    gdat.lablcurv = r'\kappa'
    gdat.lablcurvunit = ''
    gdat.lablexpo = r'\epsilon'
    gdat.lablexpounit = gdat.strgenerunit
    
    gdat.labllliktotl = r'\mathcal{L}'
    gdat.lablmedilliktotl = r'\mathcal{L}^{med}'
    
    gdat.lablasca = r'\theta_s'
    gdat.lablascaunit = gdat.lablgangunit
    gdat.lablacut = r'\theta_c'
    gdat.lablacutunit = gdat.lablgangunit
    
    gdat.lablmcut = r'M_{c,a}'
    gdat.lablmcutunit = r'$M_{\odot}$'
    
    gdat.lablspec = gdat.lablflux
    gdat.lablspecunit = gdat.lablfluxunit
    gdat.lablspecplot = gdat.lablflux
    gdat.lablspecplotunit = gdat.lablfluxunit
    gdat.lablcnts = 'C'
    gdat.lablcntsunit = ''
    gdat.labldeltllik = r'\Delta_a \ln \mathcal{L}'
    gdat.labldeltllikunit = ''
    gdat.labldiss = r'\theta_{sa}'
    gdat.labldissunit = gdat.lablgangunit
    
    gdat.labldotn = r'\Delta \theta_{pix} \hat{\theta}_a \cdot \vec{\nabla} k_d'
    gdat.labldotnunit = ''
    
    gdat.labldots = r'\vec{\alpha}_a \cdot \vec{\nabla} k_d'
    gdat.labldotsunit = ''
    
    gdat.labldotv = r'\Delta \theta_{pix} |\hat{\theta}_a \cdot \vec{\nabla} k_d|'
    gdat.labldotvunit = ''
    
    gdat.labldotm = r'\Delta \theta_{pix} \hat{\theta}_a \cdot \vec{\nabla} k_m'
    gdat.labldotmunit = ''
    
    setattr(gdat, 'labl' + gdat.namefeatsign + 'sign', getattr(gdat, 'labl' + gdat.namefeatsign))
    setattr(gdat, 'labl' + gdat.namefeatsign + 'signunit', getattr(gdat, 'labl' + gdat.namefeatsign + 'unit'))
    
    dicttemp = deepcopy(gdat.__dict__)
    for name, valu in dicttemp.iteritems():
        if name.startswith('labl') and not name.endswith('unit'):
            name = name[4:]
            labl = getattr(gdat, 'labl' + name)
            try:
                lablunit = getattr(gdat, 'labl' + name + 'unit')
                if labl.startswith('$'):
                    setattr(gdat, 'labl' + name + 'totl', '%s [%s]' % (labl, lablunit))
                else:
                    setattr(gdat, 'labl' + name + 'totl', '$%s$ [%s]' % (labl, lablunit))
            except:
                if labl.startswith('$'):
                    setattr(gdat, 'labl' + name + 'totl', '%s' % labl)
                else:
                    setattr(gdat, 'labl' + name + 'totl', '$%s$' % labl)
    
    ## legends
    if gdat.elemtype == 'lght':
        gdat.strgelem = 'PS'
    if gdat.elemtype == 'lens':
        gdat.strgelem = 'Subhalo'
    if gdat.elemtype == 'clus':
        gdat.strgelem = 'Cluster'
    
    # temp
    gdat.maxmgang = max(gdat.fittmaxmgang, gdat.truemaxmgang)

    gdat.minmspec = 1e-11
    gdat.maxmspec = 1e-7
    
    gdat.minmcnts = 0.1
    gdat.maxmcnts = 1e4
    
    gdat.minmspecplot = gdat.minmspec
    gdat.maxmspecplot = gdat.maxmspec
    
    if gdat.elemtype == 'lght':
        if gdat.exprtype == 'ferm':
            gdat.minmsind = 1.
            gdat.maxmsind = 4.
        if gdat.exprtype == 'chan':
            gdat.minmsind = -1.
            gdat.maxmsind = 3.
    
    gdat.minmexpo = 1.
    gdat.maxmexpo = 10.
    
    gdat.minmasca = 0.
    gdat.maxmasca = 0.4 / gdat.anglfact
    
    gdat.minmacut = 0.4 / gdat.anglfact
    gdat.maxmacut = 0.8 / gdat.anglfact
    
    for l in gdat.trueindxpopl:
        setattr(gdat, 'minmcmplpop%d' % l, 0.)
        setattr(gdat, 'maxmcmplpop%d' % l, 1.)
        setattr(gdat, 'scalcmplpop%d' % l, 'self')

    gdat.minmdeltllik = 1.
    gdat.maxmdeltllik = 1e3
    gdat.minmdiss = 0.
    gdat.maxmdiss = gdat.maxmgang * sqrt(2.)
    
    gdat.minmdotn = 1e0
    gdat.maxmdotn = 1e4

    gdat.minmdots = 1e0
    gdat.maxmdots = 1e4

    gdat.minmdotv = 1e4
    gdat.maxmdotv = 1e5

    gdat.minmdotm = 1e0
    gdat.maxmdotm = 1e4

    gdat.minmmcut = 5e7
    gdat.maxmmcut = 5e9
    
    gdat.minmbein = 0.
    gdat.maxmbein = 1. / gdat.anglfact
    
    # scalar variables
    gdat.minmdeflprof = 1e-3 / gdat.anglfact
    gdat.maxmdeflprof = 0.1 / gdat.anglfact

    gdat.minmfracsubh = 0.
    gdat.maxmfracsubh = 0.3
    gdat.scalfracsubh = 'self'
    
    gdat.minmmasshostbein = 1e10
    gdat.maxmmasshostbein = 1e14
    gdat.scalmasshostbein = 'self'
    
    gdat.minmmasssubhtotl = 1e8
    gdat.maxmmasssubhtotl = 1e10
    gdat.scalmasssubhtotl = 'self'

    # set up the indices of the fitting model
    retr_indxsamp(gdat)

    # construct the fitting model
    setp_fixp(gdat)
    
    for l in gdat.fittindxpopl:
        setattr(gdat, 'minmfdispop%d' % l, 0.)
        setattr(gdat, 'maxmfdispop%d' % l, 1.)
        setattr(gdat, 'scalfdispop%d' % l, 'self')
    
    for l in gdat.fittindxpopl:
        setattr(gdat, 'lablfdispop%d' % l, '$f_{%d}$' % l)

    if gdat.datatype == 'mock':
        setp_fixp(gdat, strgmodl='true')
    
    for strgmodl in gdat.liststrgmodl:
        namefixp = getattr(gdat, strgmodl + 'namefixp')
        lablfixp = getattr(gdat, strgmodl + 'lablfixp')
        lablfixpunit = getattr(gdat, strgmodl + 'lablfixpunit')
        numbfixp = getattr(gdat, strgmodl + 'numbfixp')
        scalfixp = getattr(gdat, strgmodl + 'scalfixp')
        lablfixptotl = empty(numbfixp, dtype=object)
        for k in getattr(gdat, strgmodl + 'indxfixp'):
            if lablfixpunit[k] == '':
                lablfixptotl[k] = '%s' % lablfixp[k]
            else:
                lablfixptotl[k] = '%s [%s]' % (lablfixp[k], lablfixpunit[k])
        setattr(gdat, strgmodl + 'lablfixptotl', lablfixptotl)
        for k, name in enumerate(namefixp):
            setattr(gdat, 'labl' + namefixp[k] + 'totl', lablfixptotl[k])
            setattr(gdat, 'scal' + namefixp[k], scalfixp[k])

    gdat.liststrgfeatconc = deepcopy(gdat.fittliststrgcomptotl)
    if gdat.elemtype == 'lght':
        gdat.liststrgfeatconc.remove('flux')
        gdat.liststrgfeatconc += ['spec']

    # list of scalar variable names
    gdat.listnamevarbscal = list(gdat.fittnamefixp)
    for l in gdat.trueindxpopl:
        gdat.listnamevarbscal += ['cmplpop%d' % l]
    for l in gdat.fittindxpopl:
        gdat.listnamevarbscal += ['fdispop%d' % l]
    if gdat.elemtype == 'lens':
        gdat.listnamevarbscal += ['fracsubh', 'masssubhtotl', 'masshostbein']
    gdat.numbvarbscal = len(gdat.listnamevarbscal)
    gdat.indxvarbscal = arange(gdat.numbvarbscal)
    
    # plotting factors for scalar variables
    for name in gdat.listnamevarbscal:
        if name in gdat.fittnamefixp:
            indxfixptemp = where(gdat.fittnamefixp == name)[0]
            factfixpplot = gdat.fittfactfixpplot[indxfixptemp]
            setattr(gdat, 'fact' + name + 'plot', factfixpplot)
        else:
            try:    
                getattr(gdat, 'fact' + name + 'plot')
            except:
                setattr(gdat, 'fact' + name + 'plot', 1.)

    # for each parameter in the fitting model, determine if there is a corresponding parameter in the generative model
    for k in gdat.indxvarbscal:
        try:
            temp = getattr(gdat, 'true' + gdat.namevarbscal[k])
        except:
            temp = None
        setattr(gdat, 'corr' + gdat.listnamevarbscal[k], temp)
    
	gdat.fittcorrfixp = empty(gdat.fittnumbfixp)
    for k in gdat.fittindxfixp:
        try:
            gdat.fittcorrfixp[k] = getattr(gdat, 'true' + gdat.fittnamefixp[k])
        except:
            gdat.fittcorrfixp[k] = None
 
    # copy fixp variables to individual variables
    for k, namefixp in enumerate(gdat.fittnamefixp):
        setattr(gdat, 'fittminm' + namefixp, gdat.fittminmfixp[k])
        setattr(gdat, 'fittmaxm' + namefixp, gdat.fittmaxmfixp[k])
        setattr(gdat, 'fittscal' + namefixp, gdat.fittscalfixp[k])
    
    for l in gdat.fittindxpopl:
        setattr(gdat, 'fittminmnumbpntspop%d' % l, gdat.fittminmnumbpnts[l])
        setattr(gdat, 'fittmaxmnumbpntspop%d' % l, gdat.fittmaxmnumbpnts[l])
    
    gdat.lablfeat = {}
    gdat.dictglob = {}
    gdat.lablfeatunit = {}
    gdat.lablfeattotl = {}
    setattr(gdat, 'numbbinsplot', 20)
    for strgmodl in gdat.liststrgmodl:
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        for strgfeat in liststrgfeattotl:
            # temp
            # heeeeeey
            # gdat.lablfeat[strgfeat] = getattr(gdat, 'labl' + strgfeat)
            #gdat.lablfeatunit[strgfeat] = getattr(gdat, 'labl' + strgfeat + 'unit')
            
            #if gdat.lablfeatunit[strgfeat] != '':
            #    gdat.lablfeatunit[strgfeat] = ' [%s]' % gdat.lablfeatunit[strgfeat]
            #gdat.lablfeattotl[strgfeat] = '$%s$%s' % (gdat.lablfeat[strgfeat], gdat.lablfeatunit[strgfeat])
            #if strgfeat == 'defs' or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal' or \
            #                                                                                            strgfeat == 'diss' or strgfeat == 'asca' or strgfeat == 'acut':
            #    gdat.dictglob['fact' + strgfeat + 'plot'] = gdat.anglfact
            #else:
            #    gdat.dictglob['fact' + strgfeat + 'plot'] = 1.
            
            labl = getattr(gdat, 'labl' + strgfeat)
            lablunit = getattr(gdat, 'labl' + strgfeat + 'unit')
            if lablunit != '':
                setattr(gdat, 'labl' + strgfeat + 'unit', ' [%s]' % lablunit)
            
            labltotl = '$%s$%s' % (labl, lablunit)
            setattr(gdat, 'labl' + strgfeat + 'totl', labltotl)
            if strgfeat == 'defs' or strgfeat == 'defssign' or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal' or \
                                                                                                        strgfeat == 'diss' or strgfeat == 'asca' or strgfeat == 'acut':
                setattr(gdat, 'fact' + strgfeat + 'plot', gdat.anglfact)
            else:
                setattr(gdat, 'fact' + strgfeat + 'plot', 1.)
            
            #setattr(gdat, 'numb' + strgfeat + 'plot', 20)
            
            if strgfeat == gdat.namefeatsign or strgfeat == gdat.namefeatsign + 'sign' or strgfeat == 'expo' or strgfeat == 'cnts' or \
                                                                                        strgfeat.startswith('dot') or strgfeat == 'mcut' or strgfeat == 'deltllik':
                #gdat.dictglob['scal' + strgfeat + 'plot'] = 'logt'
                setattr(gdat, 'scal' + strgfeat + 'plot', 'logt')
            else:
                #gdat.dictglob['scal' + strgfeat + 'plot'] = 'self'
                setattr(gdat, 'scal' + strgfeat + 'plot', 'self')
    
    # log-prior register
    ## indices of penalization term
    indxlpripena = 0
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    if gdat.fittmaxmnumbpntstotl > 0:
        gdat.numblpri = 1 + 7 * gdat.fittnumbpopl
    else:
        gdat.numblpri = 0

    # size of the auxiliary variable propobability density vector
    gdat.numblpau = gdat.fittmaxmnumbcomp
    
    # process index
    gdat.indxproc = arange(gdat.numbproc)

    if gdat.elemtype == 'lens':
        h5f = h5py.File(gdat.pathdata + 'inpt/adis.h5','r')
        reds = h5f['reds'][()]
        adis = h5f['adis'][()]
        adistdim = h5f['adistdim'][()]
        gdat.adisobjt = interp1d_pick(reds, adis)
        h5f.close()

        gdat.redshost = 0.2
        gdat.redssour = 1.
        gdat.adishost = gdat.adisobjt(gdat.redshost) * 1e3 # [kpc]
        gdat.adissour = gdat.adisobjt(gdat.redssour) * 1e3 # [kpc]
        gdat.adishostsour = gdat.adissour - (1. + gdat.redshost) / (1. + gdat.redssour) * gdat.adishost
        gdat.adisfact = gdat.adishost * gdat.adissour / gdat.adishostsour
        gdat.factnewtlght = 2.09e16 # Msun / kpc
        gdat.massfrombein = retr_massfrombein(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        gdat.critmden = retr_critmden(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        
        retr_axis(gdat, 'mcut', gdat.minmmcut, gdat.maxmmcut, gdat.numbbinsplot)
        retr_axis(gdat, 'bein', gdat.minmbein, gdat.maxmbein, gdat.numbbinsplot)

    # angular deviation
    # check if 1000 is too much
    # temp
    gdat.numbangl = 100
    if gdat.exprtype == 'sdyn':
        gdat.maxmangl = 1.
    if gdat.exprtype == 'ferm':
        gdat.maxmangl = 17. / gdat.anglfact
    if gdat.exprtype == 'chan':
        gdat.maxmangl = 15. / gdat.anglfact
    if gdat.exprtype == 'hubb':
        gdat.maxmangl = 2. / gdat.anglfact
    retr_axis(gdat, 'angl', 0., gdat.maxmangl, gdat.numbangl)
    retr_axis(gdat, 'anglfull', 0., 3. * gdat.maxmgang, gdat.numbangl)
    gdat.binsanglcosi = sort(cos(gdat.binsangl))
    
    # temp
    #gdat.meshbackener = meshgrid(gdat.indxback, gdat.indxener, indexing='ij')
    
    # plotting
    ## the normalized offset for text annotation of point sources in the frames
    gdat.offstext = gdat.maxmgangdata * 0.05
    ## figure size
    gdat.plotsize = 7
    ## size of the images
    gdat.sizeimag = 1.3 * gdat.plotsize
    ## text
    if gdat.datatype == 'mock':
        gdat.legdtrue = 'True'
        gdat.legdtrue = 'True'
    if gdat.datatype == 'inpt':
        gdat.legdtrue = gdat.strgcatl
        gdat.legdtrue = gdat.strgcatl

    gdat.legdtruemiss = gdat.legdtrue + ' miss'
    gdat.legdtruehits = gdat.legdtrue + ' hit'

    gdat.legdtruehost = gdat.legdtrue + ' host'
    gdat.legdtruesour = gdat.legdtrue + ' sour'
   
    # temp
    # PS indices in each population
    #gdat.indxpntspopl = []
    #for l in gdat.indxpopl:
    #    gdat.indxpntspopl.append(arange(sum(gdat.maxmnumbpnts[:l]), sum(gdat.maxmnumbpnts[:l+1])))
    
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
    retr_axis(gdat, 'oaxi', gdat.minmoaxi, gdat.maxmoaxi, gdat.numboaxi)
    gdat.binsoaxiopen = gdat.binsoaxi[:-1]
    gdat.indxoaxi = arange(gdat.numboaxi)

    # pixelization
    if gdat.datatype == 'mock':
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2

    gdat.indxchro = dict()
    gdat.indxchro['totl'] = 0
    gdat.indxchro['type'] = 1
    gdat.indxchro['prop'] = 2
    gdat.indxchro['diag'] = 3
    gdat.indxchro['save'] = 4 
    gdat.indxchro['plot'] = 5
    gdat.indxchro['proc'] = 6
    gdat.indxchro['lpri'] = 7
    gdat.indxchro['llik'] = 8
    gdat.numbchro = len(gdat.indxchro)
    
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

    if gdat.elemtype == 'lens': 
        retr_axis(gdat, 'defl', -gdat.maxmgang, gdat.maxmgang, 50)
        retr_axis(gdat, 'deflelem', -gdat.maxmgang * 1e-2, gdat.maxmgang * 1e-2, 50)

    # lensing problem setup
    ## number of deflection components to plot
    gdat.numbdeflpntsplot = 2
    gdat.numbdeflsingplot = gdat.numbdeflpntsplot + 3

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
    if gdat.elemtype == 'lens':
        gdat.adis = gdat.adishost
    
    # power spectra
    if gdat.pixltype == 'cart':
        gdat.numbwvecodim = gdat.numbsidecart
        gdat.minmanglodim = 0.
        gdat.maxmanglodim = 2. * gdat.maxmgang
        gdat.minmmpolodim = 0.
        gdat.maxmmpolodim = 1. / 2. / gdat.sizepixl
        retr_axis(gdat, 'anglodim', gdat.minmanglodim, gdat.maxmanglodim, gdat.numbsidecart, invr=True)
        retr_axis(gdat, 'mpolodim', gdat.minmmpolodim, gdat.maxmmpolodim, gdat.numbsidecart / 2)
        if gdat.elemtype == 'lens':
            gdat.minmwvecodim = gdat.minmmpolodim / gdat.adis
            gdat.maxmwvecodim = gdat.maxmmpolodim / gdat.adis
            gdat.minmwlenodim = gdat.minmanglodim * gdat.adis
            gdat.maxmwlenodim = gdat.maxmanglodim * gdat.adis
            retr_axis(gdat, 'wvecodim', gdat.minmwvecodim, gdat.maxmwvecodim, gdat.numbsidecart / 2)
            retr_axis(gdat, 'wlenodim', gdat.minmwlenodim, gdat.maxmwlenodim, gdat.numbsidecart, invr=True)
            gdat.meanwveclgal, gdat.meanwvecbgal = meshgrid(gdat.meanwvecodim, gdat.meanwvecodim, indexing='ij')
            gdat.meanwvec = sqrt(gdat.meanwveclgal**2 + gdat.meanwvecbgal**2)

    # element parameter vector indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompflux = 2
    gdat.indxcompsind = 3
    gdat.indxcompcurv = 4
    gdat.indxcompexpo = 4
    # temp
    #gdat.indxcomp = [[] for l in gdat.indxpopl]
    #for l in gdat.indxpopl:
    #    gdat.indxcomp[l] = arange(gdat.fittnumbcomp[l])
    #gdat.indxpnts = []
    #for l in gdat.indxpopl:
    #    gdat.indxpnts.append(arange(gdat.maxmnumbpnts[l]))

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
    for strgmodl in ['fitt', 'true']:
        backflux = []
        back = getattr(gdat, strgmodl + 'back')

        for c in getattr(gdat, strgmodl + 'indxback'):
            if isinstance(back[c], float):
                if gdat.pixltype == 'heal':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + back[c]
                if gdat.pixltype == 'cart':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + back[c]
            elif isinstance(back[c], ndarray) and back[c].ndim == 1:
                if gdat.pixltype == 'heal':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + back[c][:, None, None]
                if gdat.pixltype == 'cart':
                    backfluxtemp = zeros((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + back[c][:, None, None]
            else:
                path = gdat.pathinpt + back[c]
                backfluxtemp = pf.getdata(path)
                if gdat.pixltype == 'cart':
                    backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
            
            backflux.append(backfluxtemp)

        setattr(gdat, strgmodl + 'backflux', backflux)
    
    # only include desired energy and PSF class bins 
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    if gdat.correxpo:
        gdat.expo = gdat.expo[gdat.indxcubeincl]
    
    for strgmodl in ['fitt', 'true']:
        backflux = getattr(gdat, strgmodl + 'backflux')
        for c in getattr(gdat, strgmodl + 'indxback'):
            backflux[c] = backflux[c][gdat.indxcubeincl]
            if amin(backflux[c]) <= 0.:
                raise Exception('Background templates must be positive.')
        setattr(gdat, strgmodl + 'backflux', backflux)
  
    if gdat.pixltype == 'cart':
        
        for strgmodl in ['fitt', 'true']:
            backfluxcart = []
            backflux = getattr(gdat, strgmodl + 'backflux')
            for c in getattr(gdat, strgmodl + 'indxback'):
                backfluxcart.append(backflux[c].reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt)))
            setattr(gdat, strgmodl + 'backfluxcart', backfluxcart)

        gdat.expocart = gdat.expo.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))

    if gdat.datatype == 'inpt':
        gdat.exprdataflux = gdat.exprdataflux[gdat.indxcubeincl]

    # exclude voxels with vanishing exposure
    if gdat.correxpo:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.indxpixlrofi = intersect1d(gdat.indxpixlrofi, where(gdat.expo[i, :, m] > 0.)[0])
    gdat.indxcuberofi = meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.indxcube = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
    
    # plotting
    gdat.minmexpomaps = amin(gdat.expo)
    gdat.maxmexpomaps = amax(gdat.expo)

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
    
    for strgmodl in ['fitt', 'true']:
        backflux = getattr(gdat, strgmodl + 'backflux')
        for c in getattr(gdat, strgmodl + 'indxback'):
            backflux[c] = backflux[c][gdat.indxcuberofi]
        setattr(gdat, strgmodl + 'backflux', backflux)

    # temp
    #if gdat.elemtype == 'lens':
    #    gdat.backfluxlens = gdat.backflux[0][0, :, 0].reshape((gdat.numbsidecart, gdat.numbsidecart))
    
    if gdat.correxpo:
        for strgmodl in ['fitt', 'true']:
            backcnts = []
            backcntstotl = zeros_like(gdat.expo)
            for c in getattr(gdat, strgmodl + 'indxback'):
                backcntstemp = retr_cntsmaps(gdat, backflux[c])
                backcnts.append(backcntstemp)
                backcntstotl[:] += backcntstemp 
            setattr(gdat, strgmodl + 'backcnts', backcnts)
            setattr(gdat, strgmodl + 'backcntstotl', backcntstotl)
    
    if gdat.evalcirc != 'full' and gdat.numbpixl * gdat.fittmaxmnumbpntstotl < 1e5:
        gdat.calcerrr = True
    else:
        gdat.calcerrr = False
   
    if gdat.elemtype == 'lens':
        gdat.evalpsfnkern = False
    else:
        gdat.evalpsfnkern = True

    if gdat.evalpsfnkern:
        gdat.exprpsfn = retr_psfn(gdat, gdat.exprpsfp, gdat.indxener, gdat.binsangl, gdat.exprpsfntype, gdat.binsoaxi, gdat.exproaxitype)
    
    if gdat.evalcirc != 'full':
        
        if gdat.evalcirc == 'locl':
            gdat.numbprox = 3
        
        gdat.indxprox = arange(gdat.numbprox)
        gdat.binsprox = logspace(log10(gdat.fittminmflux), log10(gdat.fittmaxmflux), gdat.numbprox + 1)
        
        if gdat.evalcirc == 'locl':
            # determine the maximum angle at which the PS flux map will be computed
            gdat.maxmangleval = empty(gdat.numbprox)
            for h in gdat.indxprox:
                if gdat.specfraceval == 0:
                    gdat.maxmangleval[h] = 3. * gdat.maxmgang
                else:  
                    frac = min(1e-2, gdat.specfraceval * gdat.binsprox[0] / gdat.binsprox[h+1])
                    psfnwdth = retr_psfnwdth(gdat, gdat.exprpsfn, frac)
                    gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                    gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
        
        if gdat.exprtype == 'chan':
            gdat.maxmangleval = maximum(gdat.maxmangleval, array([7., 5., 3.]) / gdat.anglfact)
       
        if gdat.evalcirc == 'bein':
            gdat.maxmangleval = 10. * gdat.binsprox[1:]
            #gdat.maxmangleval = array([1. / gdat.anglfact]) # 4 * gdat.binsprox[1:]

        if gdat.elemtype == 'lght' and gdat.maxmangl - amax(gdat.maxmangleval) < 1.1 * sqrt(2) * (gdat.maxmgang - gdat.maxmgangdata):
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
                                                                                                            1e2 * amax(gdat.maxmangleval), gdat.numbprox)
    
    if gdat.evalcirc != 'full':

        # make a look-up table of nearby pixels for each pixel
        path = gdat.pathprox + 'indxprox_%08d_%s_%0.4g_%0.4g_%04d.p' % (gdat.numbpixl, gdat.pixltype, 1e2 * amin(gdat.maxmangleval), \
                                                                                                            1e2 * amax(gdat.maxmangleval), gdat.numbprox)
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
            gdat.indxpixlprox = [[] for h in range(gdat.numbprox)]
            cntrsave = -1.
            # temp
            for j in gdat.indxpixl:
                dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
                dist[j] = 0.
                for h in range(gdat.numbprox):
                    indxpixlproxtemp = where(dist < gdat.maxmangleval[h])[0]
                    if indxpixlproxtemp.size > 1e4:
                        indxpixlproxtemp = -1
                    gdat.indxpixlprox[h].append(indxpixlproxtemp)
                cntrsave = tdpy.util.show_prog(j, gdat.indxpixl.size, cntrsave)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()

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
    gdat.alphpnts = 1.
    gdat.alphmaps = 1.
    
    # number of colorbar ticks in the maps
    gdat.numbtickcbar = 11
    
    # turn off relevant proposal types
    gdat.indxfixpprop = []
    for k, strg in enumerate(gdat.fittnamefixp):
        if k in gdat.fittindxfixpnumbpnts:
            thisbool = False
        else:
            if k in gdat.fittindxfixpmeanpnts:
                strgtemp = 'meanpnts'
            if k in gdat.fittindxfixpdist:
                strgtemp = 'dist'
            if k in gdat.fittindxfixppsfp:
                strgtemp = 'psfp'
            if k in gdat.fittindxfixpbacp:
                strgtemp = 'bacp'
            if k in gdat.fittindxfixplenp:
                strgtemp = 'lenp'
            thisbool = getattr(gdat, 'prop' + strgtemp)
        
        if thisbool:
            gdat.indxfixpprop.append(getattr(gdat, 'fittindxfixp' + strg))
    gdat.indxfixpprop = array(gdat.indxfixpprop) 
    gdat.numbfixpprop = gdat.indxfixpprop.size
    gdat.indxprop = arange(gdat.numbfixpprop)
    
    gdat.indxstdplgal = gdat.numbfixpprop
    gdat.indxstdpbgal = gdat.numbfixpprop + 1
    if gdat.elemtype == 'lght':
        gdat.indxstdpflux = gdat.numbfixpprop + 2
        if gdat.numbener > 1:
            gdat.indxstdpsind = gdat.numbfixpprop + 3
            gdat.indxstdpcurv = gdat.numbfixpprop + 4
            gdat.indxstdpexpo = gdat.numbfixpprop + 4
    if gdat.elemtype == 'lens':
        gdat.indxstdpdefs = gdat.numbfixpprop + 2
        gdat.indxstdpasca = gdat.numbfixpprop + 3
        gdat.indxstdpacut = gdat.numbfixpprop + 4
    if gdat.elemtype == 'clus':
        gdat.indxstdpnobj = gdat.numbfixpprop + 2
    # temp
    gdat.indxstdpampl = getattr(gdat, 'indxstdp' + gdat.namecompampl)
    gdat.numbstdp = gdat.numbfixpprop + gdat.fittmaxmnumbcomp
    gdat.numbstdpfixp = gdat.numbfixpprop
    gdat.strgstdp = concatenate((array(gdat.fittlablfixp)[gdat.indxfixpprop], gdat.fittliststrgcomptotl))
    gdat.strgstdp = list(gdat.strgstdp)
    gdat.namestdp = concatenate((array(gdat.fittnamefixp)[gdat.indxfixpprop], gdat.fittliststrgcomptotl))
    gdat.indxstdp = arange(gdat.numbstdp)
    gdat.indxstdpfixp = arange(gdat.numbfixpprop)
    gdat.indxstdpcomp = setdiff1d(gdat.indxstdp, gdat.indxstdpfixp)
    gdat.indxparaprop = gdat.indxfixpprop

    # proposal scale indices for each parameter
    indxpntsfull = [range(gdat.fittmaxmnumbpnts[l]) for l in gdat.fittindxpopl]
    indxsampcomp = retr_indxsampcomp(gdat, indxpntsfull, 'fitt')
    gdat.indxstdppara = zeros(gdat.fittnumbpara, dtype=int) - 1
    cntr = 0
    gdat.indxparaprop = zeros(gdat.numbfixpprop, dtype=int)
    for k in gdat.fittindxpara:
        if k in gdat.indxfixpprop:
            gdat.indxstdppara[k] = cntr
            gdat.indxparaprop[cntr] = k
            cntr += 1
        for l in gdat.fittindxpopl:
            for strgcomp in gdat.fittliststrgcomp[l]:
                if k in indxsampcomp[strgcomp][l]:
                    gdat.indxstdppara[k] = getattr(gdat, 'indxstdp' + strgcomp)

    # for the fitting model, define proposal type indices
    dicttemptemp = deepcopy(gdat.__dict__) 
    for name, valu in dicttemptemp.iteritems():
        if name.startswith('fittindxfixp') and name != 'fittindxfixp':
            indxfixp = valu
            indxstdp = gdat.indxstdppara[indxfixp]
            setattr(gdat, 'indxstdp' + name[12:], indxstdp)
    
    # proposal scale
    if gdat.elemtype == 'lens':
        gdat.stdvstdp = 1e-2 + zeros(gdat.numbstdp)
        gdat.stdvstdp[gdat.indxstdplgalhost] = 1e-3
        gdat.stdvstdp[gdat.indxstdplgalhost] = 1e-3
        gdat.stdvstdp[gdat.indxstdpbeinhost] = 1e-3
        gdat.stdvstdp[gdat.indxstdpspechost] = 1e-3
        gdat.stdvstdp[gdat.indxstdpcomp] = 5e-2
    else:
        if gdat.exprtype == 'ferm':
            gdat.stdvstdp = 4e-5 + zeros(gdat.numbstdp)
            gdat.stdvstdp[gdat.indxstdppara[gdat.fittindxfixpmeanpnts]] = 1e-3
            gdat.stdvstdp[gdat.indxstdppara[gdat.fittindxfixpdist]] = 1e-3
            gdat.stdvstdp[gdat.indxstdpcomp] = 1e-3
            gdat.stdvstdp[gdat.indxstdpflux] = 1e-3
        if gdat.exprtype == 'chan':
            gdat.stdvstdp = 1e-2 + zeros(gdat.numbstdp)
            # temp
            gdat.stdvstdp[gdat.fittindxfixphypr+gdat.fittnumbpopl] = 1e-2
            gdat.stdvstdp[gdat.indxstdpcomp] = 1e-4
            gdat.stdvstdp[gdat.indxstdpflux] = 1e-2
        if gdat.exprtype == 'sdyn':
            gdat.stdvstdp = 1e-3 + zeros(gdat.numbstdp)

    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = retr_cntsmaps(gdat, gdat.exprdataflux)

    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        gdat.limsangl = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        gdat.limspsfn = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.trueoaxitype:
                    psfn = gdat.exprpsfn[i, :, m, 0]
                else:
                    psfn = gdat.exprpsfn[i, :, m]
                maxmpsfn = amax(psfn)
                gdat.limsangl[i][m] = [0., gdat.binsangl[amax(where(psfn > 1e-6 * maxmpsfn)[0])] * gdat.anglfact]
                gdat.limspsfn[i][m] = [maxmpsfn * 1e-6, maxmpsfn]
       
    # spatially averaged background flux 
    for strgmodl in ['fitt', 'true']:
        
        numbback = getattr(gdat, strgmodl + 'numbback')
        backflux = getattr(gdat, strgmodl + 'backflux')
        indxback = getattr(gdat, strgmodl + 'indxback')
        
        backfluxmean = zeros((numbback, gdat.numbener))
        for c in indxback:
            for i in gdat.indxener:
                if gdat.correxpo:
                    backfluxmean[c, i] += sum(backflux[c][i, :, :] * gdat.expo[i, :, :]) / sum(gdat.expo[i, :, :])
                else:
                    backfluxmean[c, i] += sum(backflux[c][i, :, 0])
        setattr(gdat, strgmodl + 'backfluxmean', backfluxmean)
    
    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(10000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)
    gdat.indxcubesave = meshgrid(gdat.indxener, gdat.indxpixlsave, gdat.indxevtt, indexing='ij')
   
    if gdat.verbtype > 1 and boolinitsetp:
        print 'fixp'
        # temp
        for strgmodl in ['fitt', 'true']:
            print 'strgmodl'
            print strgmodl
            if strgmodl == 'fitt':
                listfeat = ('name', 'strg', 'scal', 'minm', 'maxm', 'true')
                print '%20s%25s%5s%20s%20s%20s' % listfeat
            else:
                listfeat = ('name', 'strg', 'scal', 'minm', 'maxm', 'sampvarb')
                print '%20s%25s%5s%20s%20s%20s' % listfeat
            for k in getattr(gdat, strgmodl + 'indxfixp'):
                namefixp = getattr(gdat, strgmodl + 'namefixp')[k]
                lablfixp = getattr(gdat, strgmodl + 'lablfixp')[k]
                scalfixp = getattr(gdat, strgmodl + 'scalfixp')[k]
                minmfixp = getattr(gdat, strgmodl + 'minmfixp')[k]
                maxmfixp = getattr(gdat, strgmodl + 'maxmfixp')[k]
                print '%20s%25s%5s%20.6g%20.6g' % (namefixp, lablfixp, scalfixp, minmfixp, maxmfixp)
    
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        gdat.exprfwhm = 2. * retr_psfnwdth(gdat, gdat.exprpsfn, 0.5)
        gdat.stdvspatprio = amax(gdat.exprfwhm)
    if gdat.elemtype == 'lens':
        gdat.stdvspatprio = amax(gdat.exprpsfp)
    
    # proposals
    # parameters not subject to proposals
    gdat.indxfixpiact = setdiff1d(gdat.fittindxfixp, gdat.indxfixpprop)
    gdat.numbfixpiact = gdat.indxfixpiact.size
    gdat.indxiact = arange(gdat.numbfixpiact)
    
    gdat.lablproptype = array([])
    gdat.legdproptype = array([])
    gdat.nameproptype = array([])
   
    if gdat.probtran == None:
        if gdat.fittnumbtrap > 0:
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
       
    cntr = tdpy.util.cntr()
    if gdat.fittnumbtrap > 0.:
        
        if gdat.propwithsing:
            for k in gdat.indxstdp:    
                gdat.indxproptypewith = cntr.incr()
                gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{W}_{%d}$' % k)
                gdat.legdproptype = append(gdat.legdproptype, gdat.fittnamepara[where(gdat.indxstdppara == k)[0][0]])
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
    
    # data structures for JIT
    gdat.deflelem = zeros((gdat.numbpixl, 2))

    if gdat.datatype != 'mock':
        proc_datacnts(gdat)
    
   
def setpfinl(gdat, boolinitsetp=False):
    
    if 'grad' in gdat.fittspatdisttype:
        gdat.fittpdfnspatpriotemp = empty((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
        for k in range(gdat.numbsidecart + 1):
            for n in range(gdat.numbsidecart + 1):
                # temp
                temp = retr_dots(gdat, gdat.truelenscnts[0, :, :, 0], gdat.binslgalcart[k], gdat.binsbgalcart[n], 1., \
                                                                                        gdat.trueascadistmeanpop0, gdat.trueacutdistmeanpop0, gdat.indxpixl, absv=True)

                #temp /= amax(temp)
                gdat.fittpdfnspatpriotemp[k, n] = temp**gdat.dotnpowr

        gdat.fittlpdfspatprio, gdat.fittlpdfspatprioobjt = retr_spatprio(gdat, gdat.fittpdfnspatpriotemp)
        gdat.fittlpdfspatpriointp = gdat.fittlpdfspatprioobjt(gdat.bgalcart, gdat.lgalcart)
        
        print 'gdat.dotnpowr'
        print gdat.dotnpowr
        print 'gdat.fittlpdfspatpriointp'
        summgene(gdat.fittlpdfspatpriointp)
        print

    # plot settings
    ## upper limit of histograms
    if gdat.datatype == 'inpt':
        gdat.limtpntshist = [0.5, 10**ceil(log10(gdat.fittmaxmnumbpntstotl))]
    
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

    if gdat.jitt:
        gdat.llik = empty_like(gdat.datacnts)
        modlcnts = zeros_like(gdat.datacnts)
        gdat.numbthre = 4
        gdat.indxthre = arange(gdat.numbthre)
        gdat.sizechun = (gdat.numbpixl + gdat.numbthre - 1) // gdat.numbthre
        gdat.llikchun = [gdat.llik[:, k*gdat.sizechun:(k+1)*gdat.sizechun, :] for k in gdat.indxthre]
        gdat.datacntschun = [gdat.datacnts[:, k*gdat.sizechun:(k+1)*gdat.sizechun, :] for k in gdat.indxthre]
        timefunc(None, "numpy (1 thread)", retr_llik_nump, gdat, modlcnts)
        timefunc("numba (1 thread)", retr_llik_sequ, gdat.llik, gdat.datacnts, modlcnts)
        timefunc("numba (%d threads)" % gdat.numbthre, retr_llik_mult, gdat, modlcnts)


def retr_dots(gdat, maps, lgal, bgal, defs, asca, acut, indxpixltemp, absv=False):
    
    grad = dstack((gradient(maps, gdat.sizepixl, axis=0), gradient(maps, gdat.sizepixl, axis=1)))
    defl = retr_defl(gdat, lgal, bgal, defs, 0., 0., asca=asca, acut=acut, indxpixltemp=gdat.indxpixl).reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    dotstemp = sum(grad * defl, 2)
    if absv:
        dotstemp = fabs(dotstemp)
    else:
        dotstemp = -dotstemp
    
    dots = sum(dotstemp) * gdat.sizepixl
    
    return dots


def retr_llik_nump(gdat, modlcnts):
    
    llik = gdat.datacnts * log(modlcnts) - modlcnts
    
    return llik


@jit(nopython=True, nogil=True)
def retr_llik_sequ(llik, datacnts, modlcnts):
    
    for i in range(llik.shape[0]):
        for j in range(llik.shape[1]):
            for m in range(llik.shape[2]):
                llik[i, j, m] = datacnts[i, j, m] * log(modlcnts[i, j, m]) - modlcnts[i, j, m]
    

def retr_llik_mult(gdat, modlcnts):
    
    listchun = []
    for k in gdat.indxthre:
        listchun.append([gdat.llikchun[k], gdat.datacntschun[k], modlcnts[:, k*gdat.sizechun:(k+1)*gdat.sizechun, :]])
    threads = [threading.Thread(target=retr_llik_sequ, args=chun) for chun in listchun]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()


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
    gdat.maxmdatacnts = amax(sum(gdat.datacnts, 2))
    gdat.minmdatacnts = 0.01 * gdat.maxmdatacnts
    gdat.maxmresicnts = ceil(gdat.maxmdatacnts * 0.1)
    gdat.minmresicnts = -gdat.maxmresicnts

    # plotting
    liststrgcbar = ['datacnts', 'resicnts']
    for strgcbar in liststrgcbar:
        retr_ticklabl(gdat, strgcbar)


def retr_ticklabl(gdat, strgcbar):
    
    scal = getattr(gdat, 'scal' + strgcbar)

    for strglimt in ['minm', 'maxm']:
        limt = getattr(gdat, strglimt + strgcbar)
        if scal == 'asnh':
            limt = arcsinh(limt)
        if scal == 'logt':
            limt = log10(limt)
        setattr(gdat, strglimt + strgcbar, limt)

    tick = linspace(getattr(gdat, 'minm' + strgcbar), getattr(gdat, 'maxm' + strgcbar), gdat.numbtickcbar)
    labl = empty(gdat.numbtickcbar, dtype=object)
    bins = copy(tick)
    for k in range(gdat.numbtickcbar):
        if scal == 'asnh':
            bins[k] = sinh(tick[k])
        elif scal == 'logt':
            bins[k] = 10**(tick[k])
        labl[k] = '%.3g' % bins[k]

    setattr(gdat, 'bins' + strgcbar, bins)
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

    return copy(varb)


def retr_indxsamp(gdat, strgmodl='fitt'):
    
    spectype = getattr(gdat, strgmodl + 'spectype')
    spatdisttype = getattr(gdat, strgmodl + 'spatdisttype')
    
    ampldisttype = getattr(gdat, strgmodl + 'ampldisttype')
    fluxdisttype = getattr(gdat, strgmodl + 'fluxdisttype')
    defsdisttype = getattr(gdat, strgmodl + 'defsdisttype')
    
    oaxitype = getattr(gdat, strgmodl + 'oaxitype')
    psfntype = getattr(gdat, strgmodl + 'psfntype')
    maxmnumbpnts = getattr(gdat, strgmodl + 'maxmnumbpnts') 
    numbback = len(getattr(gdat, strgmodl + 'back'))
    specback = getattr(gdat, strgmodl + 'specback')
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
    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    indxpopl = arange(numbpopl, dtype=int) 


    liststrgcomp = [[] for l in indxpopl]
    listscalcomp = [[] for l in indxpopl]
    liststrgfeatprio = [[] for l in indxpopl]
    listscalfeatprio = [[] for l in indxpopl]
    liststrgfeatodim = [[] for l in indxpopl]
    liststrgfeatcorr = [[] for l in indxpopl]
    liststrgfeatsign = [[] for l in indxpopl]
    liststrgfeat = [[] for l in indxpopl]
    for l in indxpopl:
        
        liststrgcomp[l] = ['lgal', 'bgal']
        listscalcomp[l] = ['self', 'self']
        if gdat.elemtype == 'lght':
            liststrgcomp[l] += ['flux']
        if gdat.elemtype == 'lens':
            liststrgcomp[l] += ['defs']
        if gdat.elemtype == 'clus':
            liststrgcomp[l] += ['nobj']
        listscalcomp[l] += [ampldisttype]
        
        if spatdisttype[l] == 'gang':
            liststrgfeatprio[l] += ['gang', 'aang']
            listscalfeatprio[l] += ['expo', 'self']
        if spatdisttype[l] == 'unif' or spatdisttype[l] == 'gaus' or spatdisttype[l] == 'grad':
            liststrgfeatprio[l] += ['lgal', 'bgal']
            listscalfeatprio[l] += ['self', 'self']
        if spatdisttype[l] == 'disc':
            liststrgfeatprio[l] += ['lgal', 'bgal']
            listscalfeatprio[l] += ['self', 'dexp']
        if gdat.elemtype == 'lght':
            liststrgfeatprio[l] += ['flux']
            listscalfeatprio[l] += [fluxdisttype[l]]
        if gdat.elemtype == 'lens':
            liststrgfeatprio[l] += ['defs']
            listscalfeatprio[l] += [defsdisttype[l]]
        if gdat.elemtype == 'clus':
            liststrgfeatprio[l] += ['nobj']
            listscalfeatprio[l] += [defsdisttype[l]]
        
        liststrgfeatodim[l] = ['lgal', 'bgal', 'deltllik']
        liststrgfeatodim[l] += ['gang', 'aang']
        if gdat.elemtype == 'lght':
            liststrgfeatodim[l] += ['flux']
        if gdat.elemtype == 'lens':
            liststrgfeatodim[l] += ['defs']
        if gdat.elemtype == 'clus':
            liststrgfeatodim[l] += ['nobj']
        
        if strgmodl == 'true':
            liststrgfeatodim[l] += [gdat.namefeatsign + 'sign']

        if gdat.numbener > 1 and gdat.elemtype == 'lght':
            liststrgcomp[l] += ['sind']
            listscalcomp[l] += ['gaus']
            liststrgfeatprio[l] += ['sind']
            listscalfeatprio[l] += ['gaus']
            liststrgfeatodim[l] += ['sind']
            if spectype[l] == 'curv':
                liststrgcomp[l] += ['curv']
                listscalcomp[l] += ['gaus']
                liststrgfeatprio[l] += ['curv']
                listscalfeatprio[l] += ['gaus']
                if not 'curv' in liststrgfeat[l]:
                    liststrgfeat[l] += ['curv']
            if spectype[l] == 'expo':
                liststrgcomp[l] += ['expo']
                listscalcomp[l] += ['gaus']
                liststrgfeatprio[l] += ['expo']
                listscalfeatprio[l] += ['gaus']
                if not 'expo' in liststrgfeat[l]:
                    liststrgfeat[l] += ['expo']
        if gdat.elemtype == 'lens':
            liststrgfeatodim[l] += ['mcut', 'diss', 'dotn', 'dots', 'dotv', 'dotm']
            if gdat.variasca:
                liststrgcomp[l] += ['asca']
                listscalcomp[l] += ['gaus']
                liststrgfeatodim[l] += ['asca']
                liststrgfeatprio[l] += ['asca']
                listscalfeatprio[l] += ['gaus']
            if gdat.variacut:
                liststrgcomp[l] += ['acut']
                listscalcomp[l] += ['gaus']
                liststrgfeatodim[l] += ['acut']
                liststrgfeatprio[l] += ['acut']
                listscalfeatprio[l] += ['gaus']
    
        for strgfeat in liststrgfeatodim[l]:
            if strgfeat != gdat.namefeatsign + 'sign':
                liststrgfeatcorr[l].append(strgfeat)
        
        for strgfeat in liststrgfeatodim[l]:
            if strgfeat != 'lgal' and strgfeat != 'bgal':
                liststrgfeatsign[l].append(strgfeat)
        
        if gdat.elemtype == 'lght':
            liststrgfeat[l] += ['cnts', 'spec', 'specplot']
        if gdat.elemtype == 'lens':
            liststrgfeat[l] += ['deflprof']
    
        for liststrg in [liststrgcomp[l], liststrgfeatprio[l], liststrgfeatodim[l], liststrgfeatsign[l]]:
            for strgthis in liststrg:
                if not strgthis in liststrgfeat[l]:
                    liststrgfeat[l].append(strgthis)

    # list of strings across all populations
    ## element parameters
    liststrgcomptotl = retr_listconc(liststrgcomp)
    numbcomptotl = len(liststrgcomptotl)
    indxcomptotl = arange(numbcomptotl)

    ## element parameters
    liststrgfeatpriototl = retr_listconc(liststrgfeatprio)
    liststrgfeattotl = retr_listconc(liststrgfeat)
    liststrgfeatcorrtotl = retr_listconc(liststrgfeatcorr)
    
    ## one dimensional element features
    liststrgfeatodimtotl = retr_listconc(liststrgfeatodim)

    numbelemfeat = len(liststrgfeattotl)
    
    cntr = tdpy.util.cntr()
    
    liststrgfeatdefa = deepcopy(liststrgfeattotl)
    if gdat.elemtype == 'lght':
        for strgfeat in ['sind', 'curv', 'expo']:
            if not strgfeat in liststrgfeatdefa:
                liststrgfeatdefa.append(strgfeat)

    dicttemp = {}
    if maxmnumbpntstotl > 0:
        for l in indxpopl:
            dicttemp['indxfixpnumbpntspop%d' % l] = cntr.incr()
        for l in indxpopl:
            dicttemp['indxfixpmeanpntspop%d' % l] = cntr.incr()
    
        liststrgvarb = ['gangdistscal', 'bgaldistscal', 'spatdistcons', 'defsdistslop', 'fluxdistslop', 'nobjdistslop', 'sinddistmean', 'sinddiststdv', \
                                              'curvdistmean', 'curvdiststdv', 'expodistmean', 'expodiststdv']
        if gdat.elemtype == 'lens':
            if gdat.variasca:
                liststrgvarb += ['ascadistmean', 'ascadiststdv']
            if gdat.variasca:
                liststrgvarb += ['acutdistmean', 'acutdiststdv']

        # temp
        for strgvarb in liststrgvarb:
            strgtemp = 'indxfixp' + strgvarb
            temp = zeros(numbpopl, dtype=int) - 1
            dicttemp[strgtemp] = temp

        for l in range(numbpopl):
            for k, strgfeatprio in enumerate(liststrgfeatprio[l]):
                
                if listscalfeatprio[l][k] == 'expo' or listscalfeatprio[l][k] == 'dexp':
                    dicttemp['indxfixp' + strgfeatprio + 'distscalpop%d' % l] = cntr.incr()
                    dicttemp['indxfixp' + strgfeatprio + 'distscal'][l] = dicttemp['indxfixp' + strgfeatprio + 'distscalpop%d' % l]
                
                if listscalfeatprio[l][k] == 'gaus' and strgfeatprio == 'lgal':
                    dicttemp['indxfixpspatdistconspop%d' % l] = cntr.incr()
                    dicttemp['indxfixpspatdistcons'][l] = dicttemp['indxfixpspatdistconspop%d' % l]
                    
                if listscalfeatprio[l][k] == 'powr':
                    dicttemp['indxfixp' + strgfeatprio + 'distsloppop%d' % l] = cntr.incr()
                    dicttemp['indxfixp' + strgfeatprio + 'distslop'][l] = dicttemp['indxfixp' + strgfeatprio + 'distsloppop%d' % l]
                    
                if listscalfeatprio[l][k] == 'gaus':
                    dicttemp['indxfixp' + strgfeatprio + 'distmeanpop%d' % l] = cntr.incr()
                    dicttemp['indxfixp' + strgfeatprio + 'diststdvpop%d' % l] = cntr.incr()
                    dicttemp['indxfixp' + strgfeatprio + 'distmean'][l] = dicttemp['indxfixp' + strgfeatprio + 'distmeanpop%d' % l]
                    dicttemp['indxfixp' + strgfeatprio + 'diststdv'][l] = dicttemp['indxfixp' + strgfeatprio + 'diststdvpop%d' % l]
    
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

    if gdat.elemtype == 'lens':
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
    for l in indxpopl:
        if gdat.numbener > 1:
            liststrgspep[l] += ['sind']
            if spectype[l] == 'expo':
                liststrgspep[l] += ['expo']
            if spectype[l] == 'curv':
                liststrgspep[l] = ['curv']
        numbspep[l] = len(liststrgspep[l]) 

    # number of element parameters
    numbcomp = zeros(numbpopl, dtype=int)
    for l in indxpopl:
        numbcomp[l] = len(liststrgcomp[l])
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
    if numbtrap > 0:
        indxsampcompinit = indxsamptrap[0]
    else:
        indxsampcompinit = -1
    indxpara = arange(numbpara)

    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr':
            setattr(gdat, strgmodl + strg, valu)

    for strg, valu in dicttemp.iteritems():
        # temp
        if isinstance(valu, ndarray) and strg[12:16] != 'dist':
            valu = sort(valu)
        setattr(gdat, strgmodl + strg, valu)


def setp_fixp(gdat, strgmodl='fitt'):
    
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    numbtrapcuml = getattr(gdat, strgmodl + 'numbtrapcuml')
    numbtrapcumr = getattr(gdat, strgmodl + 'numbtrapcumr')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    numbfixp = getattr(gdat, strgmodl + 'numbfixp')
    numbback = getattr(gdat, strgmodl + 'numbback')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    numbpara = getattr(gdat, strgmodl + 'numbpara')
    nameback = getattr(gdat, strgmodl + 'nameback')

    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    psfntype = getattr(gdat, strgmodl + 'psfntype')
    numbpsfptotl = getattr(gdat, strgmodl + 'numbpsfptotl')
    specback = getattr(gdat, strgmodl + 'specback')
    indxbackbacp = getattr(gdat, strgmodl + 'indxbackbacp')
    indxenerbacp = getattr(gdat, strgmodl + 'indxenerbacp')

    liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')

    # construct the fixed parameter structure
    namefixp = zeros(numbfixp, dtype=object)
    lablfixp = zeros(numbfixp, dtype=object)
    lablfixpunit = zeros(numbfixp, dtype=object)
    lablfixpunit[:] = ''
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
        
        strgvarb = strg[4:]

        if not isinstance(k, int) or not strgvarb.startswith('indxfixp') or strgvarb == 'indxfixppsfpinit' or strg[:4] != strgmodl:
            continue

        strgvarb = strgvarb[8:]
        
        namefixp[k] = strgvarb

        if strg[:-1].endswith('pop'):
            
            if numbpopl == 1:
                strgpopl = ''
                strgpoplcomm = ''
            else:
                strgpopl = '%s' % strg[-1]
                strgpoplcomm = ',%s' % strg[-1]
            
            if namefixp[k].startswith('numbpnts'):
                lablfixp[k] = '$N_{%s}$' % strgpopl
                scalfixp[k] = 'pois'
                
            if namefixp[k].startswith('meanpnts'):
                lablfixp[k] = r'$\mu_{%s}$' % strgpopl
                scalfixp[k] = 'logt'
    
            if namefixp[k].startswith('gangdistscalp'):
                lablfixp[k] = r'$\gamma_{\theta%s}$' % strgpoplcomm
                scalfixp[k] = 'self'
                lablfixpunit[k] = gdat.lablgangunit
            
            if namefixp[k].startswith('spatdistconsp'):
                lablfixp[k] = r'$\gamma_{c%s}$' % strgpoplcomm
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('bgaldistscalp'):
                lablfixp[k] = r'$\gamma_{b%s}$' % strgpoplcomm
                scalfixp[k] = 'self'
                lablfixpunit[k] = gdat.lablgangunit
            
            if namefixp[k][4:].startswith('distslopp'):
                lablfixp[k] = r'$\beta_{%s}$' % strgpopl
                if gdat.elemtype == 'lens':
                    scalfixp[k] = 'gaus'
                else:
                    scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('sinddistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablsind)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('sinddiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablsind)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('curvdistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablcurv)
                scalfixp[k] = 'atan'
            
            if namefixp[k].startswith('curvdiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablcurv)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('expodistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablexpo)
                scalfixp[k] = 'logt'
                lablfixpunit[k] = gdat.lablenerunit
            
            if namefixp[k].startswith('expodiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablexpo)
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
                n = k - getattr(gdat, strgmodl + 'indxfixpsigcene0evt0')
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
                lablfixpunit[k] = gdat.lablgangunit
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
                indxenertemp = gdat.indxenerincl[((k - getattr(gdat, strgmodl + 'indxfixpsigcene0evt0')) % (gdat.numbener * numbpsfptotl)) // numbpsfptotl]
                strgenertemp = '%s' % indxenertemp
            else:
                strgenertemp = ''
            if gdat.numbevtt > 1:
                indxevtttemp = gdat.indxevttincl[(k - getattr(gdat, strgmodl + 'indxfixpsigcene0evt0')) // (gdat.numbener * numbpsfptotl)]
                strgevtttemp = '%s' % indxevtttemp
            else:
                strgevtttemp = ''
            lablfixp[k] = r'$%s^{%s}_{%s%s}$' % (strgvarbtemp, strgcomptemp, strgenertemp, strgevtttemp)
        
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
            lablfixp[k] = '$A_{%s%s}$' % (strgbacktemp, strgenertemp)
            lablfixpunit[k] = gdat.lablfluxsoldunit
            scalfixp[k] = 'logt'
        
        if gdat.elemtype == 'lens':
            if k in getattr(gdat, strgmodl + 'indxfixplenp'):
                if strgvarb == 'lgalsour':
                    lablfixp[k] = '$l_{src}$'
                    scalfixp[k] = 'gaus'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'bgalsour':
                    lablfixp[k] = '$b_{src}$'
                    scalfixp[k] = 'gaus'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('specsour'):
                    if gdat.numbener > 1:
                        lablfixp[k] = '$f_{src,%s}$' % strg[-1]
                    else:
                        lablfixp[k] = '$f_{src}$'
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if strgvarb == 'sizesour':
                    lablfixp[k] = '$R_{e,src}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'ellpsour':
                    lablfixp[k] = r'$\epsilon_{src}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'anglsour':
                    lablfixp[k] = r'$\phi_{src}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'lgalhost':
                    lablfixp[k] = '$l_{hst}$'
                    scalfixp[k] = 'gaus'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'bgalhost':
                    lablfixp[k] = '$b_{hst}$'
                    scalfixp[k] = 'gaus'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('spechost'):
                    if gdat.numbener > 1:
                        lablfixp[k] = '$f_{hst,%s}$' % strg[-1]
                    else:
                        lablfixp[k] = '$f_{hst}$'
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if strgvarb == 'sizehost':
                    lablfixp[k] = '$R_{e,hst}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'beinhost':
                    lablfixp[k] = r'$\theta_{E,hst}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'ellphost':
                    lablfixp[k] = r'$\epsilon_{hst}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'anglhost':
                    lablfixp[k] = r'$\phi_{hst}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'sherhost':
                    lablfixp[k] = r'$\gamma_{ext}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'sanghost':
                    lablfixp[k] = r'$\phi_{ext}$'
                    scalfixp[k] = 'self'
        
        if scalfixp[k] == 'pois' or scalfixp[k] == 'self' or scalfixp[k] == 'logt' or scalfixp[k] == 'atan':
            if strgvarb.startswith('numbpnts'):
                minmfixp[k] = getattr(gdat, strgmodl + 'minm' + strgvarb[:-4])[k]
                maxmfixp[k] = getattr(gdat, strgmodl + 'maxm' + strgvarb[:-4])[k]
            else:
                minmfixp[k] = getattr(gdat, strgmodl + 'minm' + strgvarb)
                maxmfixp[k] = getattr(gdat, strgmodl + 'maxm' + strgvarb)
        if scalfixp[k] == 'gaus' or scalfixp[k] == 'eerr':
            meanfixp[k] = getattr(gdat, strgmodl + 'mean' + strgvarb)
            stdvfixp[k] = getattr(gdat, strgmodl + 'stdv' + strgvarb)
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
   
    # background templates
    listlablcompfrac = deepcopy(nameback)
    if gdat.elemtype == 'lght':
        listlablcompfrac.append('PS')
    if gdat.elemtype == 'lens':
        listlablcompfrac.append('Source')
        listlablcompfrac.append('Host')
    if gdat.elemtype == 'clus':
        listlablcompfrac.append('Uniform')
    listlablcompfracspec = deepcopy(listlablcompfrac)
    listlablcompfracspec += ['Data']
    if len(listlablcompfrac) > 1:
        listlablcompfracspec.append('Total Model')
    numblablcompfrac = len(listlablcompfrac)
    numblablcompfracspec = len(listlablcompfracspec)

    namepara = zeros(int(numbpara), dtype=object)
    namepara = zeros(numbpara, dtype=object)
    scalpara = zeros(numbpara, dtype=object)
    
    if strgmodl == 'fitt':
        for k in indxfixp:
            setattr(gdat, 'labl' + namefixp[k], lablfixp[k])
            setattr(gdat, 'fact' + namefixp[k] + 'plot', factfixpplot[k])
        
    namepara[indxfixp] = namefixp
    scalpara[indxfixp] = scalfixp
    for k in range(numbtrap):
        indxpopltemp = argmin(where(k // numbtrapcumr == 0))
        indxcomptemp = (k - numbtrapcuml[indxpopltemp]) % numbcomp[indxpopltemp]
        namepara[numbfixp+k] = liststrgcomp[indxpopltemp][indxcomptemp]
        scalpara[numbfixp+k] = listscalcomp[indxpopltemp][indxcomptemp]
    
    #listscalparatotl = list(unique(scalpara))
    #listnameparatotl = list(unique(namepara))
    # limits for the Gaussian distributed parameters
    numbstdvgaus = 2.
    #for k, namepara in enumerate(listnameparatotl):
    for k, name in enumerate(namepara):

        if scalpara[k] == 'gaus':
            if name in liststrgcomptotl:
                minm = 1e100
                maxm = -1e100
                for l in indxpopl:
                    indx = liststrgcomp[l].index(name)
                    if listscalcomp[l][indx] == 'gaus':
                        minmtemp = getattr(gdat, strgmodl + 'minm' + name + 'distmeanpop%d' % l) - numbstdvgaus * getattr(gdat, strgmodl + 'maxm' + name + 'diststdvpop%d' % l)
                        maxmtemp = getattr(gdat, strgmodl + 'maxm' + name + 'distmeanpop%d' % l) + numbstdvgaus * getattr(gdat, strgmodl + 'maxm' + name + 'diststdvpop%d' % l)
                        minm = min(minmtemp, minm)
                        maxm = max(maxmtemp, maxm)
            else:
                minm = getattr(gdat, strgmodl + 'mean' + name) - numbstdvgaus * getattr(gdat, strgmodl + 'stdv' + name)
                maxm = getattr(gdat, strgmodl + 'mean' + name) + numbstdvgaus * getattr(gdat, strgmodl + 'stdv' + name)

            setattr(gdat, strgmodl + 'minm' + name, minm)
            setattr(gdat, strgmodl + 'maxm' + name, maxm)
    
    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr':
            setattr(gdat, strgmodl + strg, valu)


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
            truestrgvarb = getattr(gdat, 'true' + strgvarb)
            if truestrgvarb == None:
                raise
        except:
            setattr(gdat, 'true' + strgvarb, valu)
 

def setp_truedefa(gdat, strgvarb, listvalu, typelimt='minmmaxm', popl=False, ener=False, evtt=False, back=False):
    
    liststrgvarb = setp_truetemp(gdat, strgvarb, popl, ener, evtt, back)

    for strgvarb in liststrgvarb:

        if typelimt == 'minmmaxm':
            try:
                minmpara = getattr(gdat, 'trueminm' + strgvarb)
                if minmpara == None:
                    raise
            except:
                setattr(gdat, 'trueminm' + strgvarb, listvalu[0])
            try:
                maxmpara = getattr(gdat, 'truemaxm' + strgvarb)
                if maxmpara == None:
                    raise
            except:
                setattr(gdat, 'truemaxm' + strgvarb, listvalu[1])
        else:
            try:
                meanpara = getattr(gdat, 'truemean' + strgvarb)
                if meanpara == None:
                    raise
            except:
                setattr(gdat, 'truemean' + strgvarb, listvalu[0])
            try:
                stdvpara = getattr(gdat, 'truestdv' + strgvarb)
                if stdvpara == None:
                    raise
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

    return fluxbrgt, fluxbrgtassc


def retr_indxoaxipnts(gdat, lgal, bgal):

    oaxi = retr_angldist(gdat, lgal, bgal, 0., 0.)
    # temp -- oaxi may need a [0]
    indxoaxipnts = digitize(oaxi, gdat.binsoaxiopen) - 1
   
    return indxoaxipnts


def init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=None, indxevttplot=None, indxpoplplot=-1):

    if strg == 'this' or gdatmodi != None:
        pathfold = gdat.pathplot + gdat.namesampdist + '/fram/'
    elif strg == 'true' or strg == '':
        pathfold = gdat.pathinit
    elif strg == 'post':
        pathfold = gdat.pathplot + gdat.namesampdist + '/finl/'
    
    figr, axis = plt.subplots(figsize=(gdat.sizeimag, gdat.sizeimag))
    
    if indxenerplot == None:
        strgener = ''
    else:
        strgener = 'ene%d' % gdat.indxenerincl[indxenerplot]
    
    if indxevttplot == None:
        strgevtt = ''
    elif indxevttplot == -1:
        strgevtt = 'evtA'
    else:
        strgevtt = 'evt%d' % gdat.indxevttincl[indxevttplot]
    
    if indxpoplplot == -1:
        strgpopl = 'popA'
    else:
        strgpopl = 'pop%d' % indxpoplplot

    if gdatmodi == None:
        strgswep = ''
    else:
        strgswep = '_swep%09d' % gdatmodi.cntrswep
    
    path = '%s%s%s%s%s%s.pdf' % (pathfold, strgplot, strgener, strgevtt, strgpopl, strgswep)
   
    axis.set_xlabel(gdat.labllgaltotl)
    axis.set_ylabel(gdat.lablbgaltotl)
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


def retr_imag(gdat, axis, maps, strg, strgcbar, thisindxener=None, thisindxevtt=-1, tdim=False):
    
    vmin = getattr(gdat, 'minm' + strgcbar)
    vmax = getattr(gdat, 'maxm' + strgcbar)
    cmap = getattr(gdat, 'cmap' + strgcbar) 
    scal = getattr(gdat, 'scal' + strgcbar) 

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
    if scal == 'logt':
        maps = log10(maps)
    
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
                                                                                                label=gdat.legdtruehits, marker='x', linewidth=2, color='g')
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, facecolor='none', \
                                                                                                    label=gdat.legdtruemiss, marker='s', linewidth=2, color='g')
    if gdat.elemtype == 'lens':
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='Model Source', marker='<', linewidth=2, color='b')
    
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='Model Host', marker='s', linewidth=2, color='b')
        if gdat.trueinfo:
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='%s Source' % gdat.legdtrue, marker='>', linewidth=2, color='g')
        
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphpnts, \
                                                                                                label='%s Host' % gdat.legdtrue, marker='D', linewidth=2, color='g')
    

    temphand, temp = axis.get_legend_handles_labels()
    numblabl = len(temp)
    
    if numblabl == 4:
        numbcols = 2
    else:
        numbcols = 3
    axis.legend(bbox_to_anchor=[0.5, 1.15], loc='center', ncol=numbcols)
        

def supr_fram(gdat, gdatmodi, strg, axis, indxpoplplot=-1):

    # true catalog
    if gdat.trueinfo:
        if indxpoplplot == -1:
            listindxpoplplot = gdat.trueindxpopl
        else:
            listindxpoplplot = [indxpoplplot]
        for l in listindxpoplplot:
        
            if gdat.truenumbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, getattr(gdat, 'true' + gdat.namefeatsign)[l][0, :])
                lgal = copy(gdat.truelgal[l][0, :])
                bgal = copy(gdat.truebgal[l][0, :])
                numbpnts = int(gdat.truenumbpnts[l])
                
                if gdatmodi != None:
                   
                    ## associations
                    ### missed
                    indx = gdatmodi.trueindxpntsasscmiss[l]
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.legdtruemiss, facecolor='none', \
                                                                                                                                marker='s', linewidth=2, color='g')
                    
                    ### hit
                    indx = gdatmodi.trueindxpntsasschits[l]
                    
                else:
                    indx = arange(lgal.size)
                
                axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, \
                                                                                        label=gdat.legdtruehits, marker='x', linewidth=2, color='g')
            
            if gdat.elemtype == 'lens':
               
                ## host
                axis.scatter(gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalhost], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalhost], \
                                                                            alpha=gdat.alphpnts, label=gdat.legdtruehits, s=300, marker='D', linewidth=2, color='g')
                axis.add_patch(plt.Circle((gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalhost], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalhost]), \
                                                                                gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbeinhost], edgecolor='g', facecolor='none', lw=2))
                
                ## source
                axis.scatter(gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalsour], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalsour], \
                                                                                    alpha=gdat.alphpnts, label=gdat.legdtruehits, s=300, marker='>', linewidth=2, color='g')
            
    # model catalog
    if indxpoplplot == -1:
        listindxpoplplot = gdat.fittindxpopl
    else:
        listindxpoplplot = [indxpoplplot]
    for l in listindxpoplplot:
        if gdatmodi != None:
            if gdat.fittnumbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namefeatsign][l]])
                lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][l]]
                bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][l]]
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphpnts, label='Sample', marker='+', linewidth=2, color='b')

            if gdat.elemtype == 'lens':
                
                ## source
                lgalsour = gdatmodi.thissampvarb[gdat.fittindxfixplgalsour]
                bgalsour = gdatmodi.thissampvarb[gdat.fittindxfixpbgalsour]
                axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, alpha=gdat.alphpnts, label='Model Source', s=300, marker='<', linewidth=2, color='b')
    
                ## host
                lgalhost = gdatmodi.thissampvarb[gdat.fittindxfixplgalhost]
                bgalhost = gdatmodi.thissampvarb[gdat.fittindxfixpbgalhost]
                beinhost = gdatmodi.thissampvarb[gdat.fittindxfixpbeinhost]
                axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, alpha=gdat.alphpnts, label='Model Host', s=300, marker='s', linewidth=2, color='b')
                axis.add_patch(plt.Circle((gdat.anglfact * lgalhost, gdat.anglfact * bgalhost), gdat.anglfact * beinhost, edgecolor='b', facecolor='none', lw=2, ls='--'))
                
    # temp
    if strg == 'post' and gdat.condcatl:
        lgal = zeros(gdat.numbprvlhigh)
        bgal = zeros(gdat.numbprvlhigh)
        sign = zeros(gdat.numbprvlhigh)
        cntr = 0
        for r in gdat.indxstkscond:
            if r in gdat.indxprvlhigh:
                lgal[cntr] = gdat.dictglob['poststkscond'][r]['lgal'][0]
                bgal[cntr] = gdat.dictglob['poststkscond'][r]['bgal'][0]
                sign[cntr] = gdat.dictglob['poststkscond'][r][gdat.namefeatsign][0]
                cntr += 1
        mrkrsize = retr_mrkrsize(gdat, sign)
        axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, label='Condensed', marker='+', linewidth=2, color='black')
        if False:
            for r in gdat.indxstkscond:
                lgal = array([gdat.dictglob['liststkscond'][r]['lgal']])
                bgal = array([gdat.dictglob['liststkscond'][r]['bgal']])
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, marker='+', linewidth=2, color='black', alpha=0.1)


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
    
    # temp -- heal does not work when the dimension of lgalfrst is 1
    if gdat.pixltype == 'heal':
        dir1 = array([lgalfrst, bgalfrst])
        dir2 = array([lgalseco, bgalseco])
        angldist = hp.rotator.angdist(dir1, dir2)
        if False:
            print 'lgalfrst'
            print lgalfrst
            print 'bgalfrst'
            print bgalfrst
            print 'lgalseco'
            print lgalseco
            print 'bgalseco'
            print bgalseco
            print 'dir1'
            print dir1
            print 'dir2'
            print dir2
            print 'angldist'
            print angldist
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
        if attr.startswith('list'):
            attrtemp = attr[4:]

            # temp
            if attrtemp == 'meanpntspop0':
                attrtemp = 'meanpnts'
            if attrtemp == 'fluxdistsloppop0':
                attrtemp = 'fluxdistslop'

            if attrtemp in gdat.fittliststrgfeatodimtotl:
                # temp
                if attrtemp == 'deltllik':
                    continue

                for l in gdat.fittindxpopl:
                    attrprim = attrtemp + 'pop%d' % l
                    listtemp = []
                    arrytemp = zeros((gdat.numbsamptotl, gdat.fittmaxmnumbpnts[l])) + nan
                    for n in range(len(valu)):
                        arrytemp[n, :len(valu[n][l])] = valu[n][l]
                    thisfile.create_dataset(attrprim, data=arrytemp)
            elif attrtemp == 'meanpnts' or attrtemp == 'fluxdistslop' or attrtemp == 'psfp' or attrtemp == 'bacp':
                thisfile.create_dataset(attrtemp, data=valu)
                
    thisfile.close()
    

def retr_deflcutf(angl, defs, asca, acut, asym=False):

    fracanglasca = angl / asca
    
    fact = ones_like(fracanglasca)
    indxlowr = where(fracanglasca < 1.)[0]
    indxuppr = where(fracanglasca > 1.)[0]
    fact[indxlowr] = arccosh(1. / fracanglasca[indxlowr]) / sqrt(1. - fracanglasca[indxlowr]**2)
    fact[indxuppr] = arccos(1. / fracanglasca[indxuppr]) / sqrt(fracanglasca[indxuppr]**2 - 1.)
    
    deflcutf = defs / fracanglasca / (1. - log(2.))
    if asym:
        deflcutf *= fact + log(fracanglasca / 2.)
    else:
        fracacutasca = acut / asca
        factcutf = fracacutasca**2 / (fracacutasca**2 + 1)**2 * ((fracacutasca**2 + 1. + 2. * (fracanglasca**2 - 1.)) * fact + \
                pi * fracacutasca + (fracacutasca**2 - 1.) * log(fracacutasca) + sqrt(fracanglasca**2 + fracacutasca**2) * (-pi + (fracacutasca**2 - 1.) / fracacutasca * \
                log(fracanglasca / (sqrt(fracanglasca**2 + fracacutasca**2) + fracacutasca))))
        deflcutf *= factcutf
       
    return deflcutf


def initchro(gdat, gdatmodi, name):

    if gdatmodi != None:
        gdatmodi.thischro[gdat.indxchro[name]] = gdat.functime()
    

def stopchro(gdat, gdatmodi, name):
    
    if gdatmodi != None:
        gdatmodi.thischro[gdat.indxchro[name]] = gdat.functime() - gdatmodi.thischro[gdat.indxchro[name]]


def retr_defl(gdat, *listargs, **dictargskeyw):
    
    defl = retr_defl_jitt(gdat.indxpixl, gdat.lgalgrid, gdat.bgalgrid, *listargs, **dictargskeyw)

    return defl

   
#@jit(nopython=True, nogil=True)
def retr_defl_jitt(indxpixl, lgalgrid, bgalgrid, lgal, bgal, bein, ellp, angl, rcor=0., asca=None, acut=None, indxpixltemp=None):
    
    if indxpixltemp == None:
        indxpixltemp = indxpixl
    
    # translate the grid
    lgaltran = lgalgrid[indxpixltemp] - lgal
    bgaltran = bgalgrid[indxpixltemp] - bgal
    
    if acut != None:
        angl = sqrt(lgaltran**2 + bgaltran**2)
        defl = retr_deflcutf(angl, bein, asca, acut)
        
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
   
    defl = vstack((defllgal, deflbgal)).T
    
    return defl


def retr_lpripowrdist(gdat, gdatmodi, strgmodl, comp, strgcomp, sampvarb, l):
    
    distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]

    minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
    maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
    lpri = sum(log(pdfn_powr(comp, minm, maxm, distslop)))

    return lpri


def retr_lpriigamdist(gdat, gdatmodi, strgmodl, comp, strgcomp, sampvarb, l):
    
    distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]

    cutf = getattr(gdat, strgmodl + 'cutf' + strgcomp)
    lpri = sum(log(igam_powr(comp, distslop, cutf)))

    return lpri


def retr_lprigausdist(gdat, gdatmodi, strgmodl, comp, strgcomp, sampvarb, l):
    
    distmean = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distmean')[l]]
    diststdv = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'diststdv')[l]]
    lpri = sum(log(pdfn_gaus(comp, distmean, diststdv)))
    
    return lpri


def traptdim(gdat, arry):
    
    s1 = arry[0, 0] + arry[-1, 0] + arry[0, -1] + arry[-1, -1]
    s2 = sum(arry[1:-1, 0]) + sum(arry[1:-1, -1]) + sum(arry[0, 1:-1]) + sum(arry[-1, 1:-1])
    s3 = sum(arry[1:-1, 1:-1])
    summ = (s1 + 2*s2 + 4*s3) * gdat.apixmodl
    
    return summ


def retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons=None):
    
    pdfnspatprio = pdfnspatpriotemp
    if spatdistcons != None:
        pdfnspatprio += spatdistcons

    summ = traptdim(gdat, pdfnspatprio)
    pdfnspatprio /= summ
    print 'pdfnspatprio'
    summgene(pdfnspatprio)
    lpdfspatprio = log(pdfnspatprio)
    lpdfspatprioobjt = sp.interpolate.RectBivariateSpline(gdat.binsbgalcart, gdat.binslgalcart, lpdfspatprio)
    
    return lpdfspatprio, lpdfspatprioobjt


@jit(nopython=True, nogil=True)
def retr_deflelem_jitt(deflelem, indxpixl, lgalgrid, bgalgrid, numbpntsconc, lgalconc, bgalconc, defsconc, ascaconc, acutconc):
    
    for k in range(numbpntsconc):
        deflelem[:] += retr_defl_jitt(indxpixl, lgalgrid, bgalgrid, lgalconc[k], bgalconc[k], defsconc[k], 0., 0., asca=ascaconc[k], acut=acutconc[k])


def proc_samp(gdat, gdatmodi, strg, raww=False, fast=False, lprionly=False):

    initchro(gdat, gdatmodi, 'proc')

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
        strgmodl = strg
    else:
        strgmodl = 'fitt'
    
    # temp
    indxfixpmeanpnts = getattr(gdat, strgmodl + 'indxfixpmeanpnts')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    spatdisttype = getattr(gdat, strgmodl + 'spatdisttype')
    ampldisttype = getattr(gdat, strgmodl + 'ampldisttype')
    if gdat.elemtype == 'lght':
        fluxdisttype = getattr(gdat, strgmodl + 'fluxdisttype')
        minmflux = getattr(gdat, strgmodl + 'minmflux')
    if gdat.elemtype == 'lens':
        defsdisttype = getattr(gdat, strgmodl + 'defsdisttype')
        minmdefs = getattr(gdat, strgmodl + 'minmdefs')
    if gdat.elemtype == 'clus':
        nobjdisttype = getattr(gdat, strgmodl + 'nobjdisttype')
        minmnobj = getattr(gdat, strgmodl + 'minmnobj')
        binsnobj = getattr(gdat, strgmodl + 'binsnobj')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    numbpopl = getattr(gdat, strgmodl + 'numbpopl')
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    spectype = getattr(gdat, strgmodl + 'spectype')
    liststrgfeatdefa = getattr(gdat, strgmodl + 'liststrgfeatdefa')
    liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
    liststrgfeatodimtotl = getattr(gdat, strgmodl + 'liststrgfeatodimtotl')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
    if gdat.evalpsfnkern:
        oaxitype = getattr(gdat, strgmodl + 'oaxitype')
    
    # common dictionary
    dicttemp = {}
           
    # grab the sample vector
    sampvarb = getattr(gdatobjt, strg + 'sampvarb')
    
    psfp = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfp')]
    bacp = sampvarb[getattr(gdat, strgmodl + 'indxfixpbacp')]
    
    if numbtrap > 0:
        indxpntsfull = list(getattr(gdatobjt, strg + 'indxpntsfull'))
        # temp -- this may slow down execution
        indxsampcomp = retr_indxsampcomp(gdat, indxpntsfull, strgmodl)
        setattr(gdatobjt, strg + 'indxsampcomp', indxsampcomp)
        indxfixpnumbpnts = getattr(gdat, strgmodl + 'indxfixpnumbpnts')
    
        numbpnts = sampvarb[indxfixpnumbpnts].astype(int)
        
        for strgfeat in liststrgfeatdefa:
            dicttemp[strgfeat] = [[] for l in range(numbpopl)]
        for l in indxpopl:
            for strgcomp in liststrgcomp[l]:
                dicttemp[strgcomp][l] = sampvarb[indxsampcomp[strgcomp][l]]
        # temp
        if gdat.elemtype == 'lens':
            for l in gdat.fittindxpopl:
                if gdat.variasca:
                    indx = where(sampvarb[indxsampcomp['acut'][l]] < 0.)[0]
                    if indx.size > 0:
                        print 'Acut went negative'
                        sampvarb[indxsampcomp['acut'][l]][indx] = 1e-3 * gdat.anglfact
                if gdat.variacut:
                    indx = where(sampvarb[indxsampcomp['asca'][l]] < 0.)[0]
                    if indx.size > 0:
                        print 'Asca went negative'
                        sampvarb[indxsampcomp['asca'][l]][indx] = 1e-3 * gdat.anglfact

        if gdat.elemtype == 'lght':
            for l in range(numbpopl):
                dicttemp['spec'][l] = retr_spec(gdat, dicttemp['flux'][l], dicttemp['sind'][l], dicttemp['curv'][l], dicttemp['expo'][l], spectype=spectype[l])
        
        for strgfeat in gdat.liststrgfeatconc:
            if strgfeat == 'spec':
                dicttemp['specconc'] = concatenate(dicttemp['spec'], axis=1)
            else:
                dicttemp[strgfeat + 'conc'] = concatenate(dicttemp[strgfeat])
        
        numbpntsconc = dicttemp['lgalconc'].size
    
    if strg == 'next' and gdat.verbtype > 1:
        show_samp(gdat, gdatmodi)
    
    # log-prior
    initchro(gdat, gdatmodi, 'lpri')

    lpri = zeros(gdat.numblpri)
    if numbtrap > 0:
    
        if 'gaus' in spatdisttype or 'grad' in spatdisttype:

            if 'gaus' in spatdisttype:
                pdfnspatpriotemp = getattr(gdat, strgmodl + 'pdfnspatpriotemp')
                spatdistcons = sampvarb[getattr(gdat, strgmodl + 'indxfixpspatdistcons')]
                lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons)
                lpdfspatpriointp = lpdfspatprioobjt(gdat.bgalcart, gdat.lgalcart)
                
                # temp
                lpdfspatpriointp = lpdfspatpriointp.T
                
                setattr(gdatobjt, strg + 'lpdfspatpriointp', lpdfspatpriointp)
                setattr(gdatobjt, strg + 'lpdfspatprioobjt', lpdfspatprioobjt)
        
            else:
                lpdfspatprioobjt = gdat.fittlpdfspatprioobjt
                
        meanpnts = sampvarb[indxfixpmeanpnts]
        
        for l in gdat.fittindxpopl:
            lpri[0] -= 0.5 * gdat.priofactdoff * numbcomp[l] * numbpnts[l]
            lpri[1+0*numbpopl+l] = retr_probpois(numbpnts[l], meanpnts[l])
            for k, strgcomp in enumerate(liststrgcomp[l]):
                if strgcomp == 'lgal':
                    if spatdisttype[l] == 'unif':
                        lpri[1+numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang)
                        lpri[1+2*numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang)
                    elif spatdisttype[l] == 'disc':
                        lpri[1+numbpopl+l] = -numbpnts[l] * log(2. * gdat.maxmgang)
                        indxfixpbgaldistscal = getattr(gdat, strgmodl + 'indxfixpbgaldistscalpop%d' % l)
                        lpri[1+2*numbpopl+l] = sum(log(pdfn_dexp(dicttemp['bgal'][l], gdat.maxmgang, sampvarb[indxfixpbgaldistscal]))) 
                    elif spatdisttype[l] == 'gang':
                        gang = retr_gang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                        indxfixpgangdistscal = getattr(gdat, strgmodl + 'indxfixpgangdistscalpop%d' % l)
                        lpri[1+numbpopl+l] = sum(log(pdfn_expo(gang, gdat.maxmgang, sampvarb[indxfixpgangdistscal]))) 
                        lpri[1+2*numbpopl+l] = -numbpnts[l] * log(2. * pi) 
                    elif spatdisttype[l] == 'gaus' or spatdisttype[l] == 'grad':
                        lpri[1+numbpopl+l] = sum(lpdfspatprioobjt(dicttemp['bgal'][l], dicttemp['lgal'][l], grid=False))
                    
                elif strgcomp == 'bgal':
                    continue 
                else: 
                    if listscalcomp[l][k] == 'powr':
                        lpri[1+(k+1)*numbpopl+l] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, dicttemp[strgcomp][l], strgcomp, sampvarb, l)
                    if listscalcomp[l][k] == 'gaus':
                        lpri[1+(k+1)*numbpopl+l] = retr_lprigausdist(gdat, gdatmodi, strgmodl, dicttemp[strgcomp][l], strgcomp, sampvarb, l)
            
        if strg == 'this':
            gdatmodi.thislpripena = lpri[0]
        
        if strg == 'next' and gdatmodi.proptran:
            
            gdatmodi.thislpau = zeros(gdat.numblpau)
            
            if gdatmodi.propbrth or gdatmodi.propsplt:
                sampvarbtemp = gdatmodi.nextsampvarb
            if gdatmodi.propdeth or gdatmodi.propmerg:
                sampvarbtemp = gdatmodi.thissampvarb
        
            if gdatmodi.propbrth or gdatmodi.propdeth:
            
                for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi]):
                    if strgcomp == 'lgal' or strgcomp == 'bgal':
                        gdatmodi.thislpau[gdat.fittmaxmnumbcomp*l+k] = -log(2. * gdat.fittmaxmgang)
                    if listscalcomp[gdatmodi.indxpoplmodi][k] == 'powr':
                        gdatmodi.thislpau[k] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                           strgcomp, gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                    if listscalcomp[gdatmodi.indxpoplmodi][k] == 'gaus':
                        gdatmodi.thislpau[k] = retr_lprigausdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                           strgcomp, gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                
            if gdatmodi.propbrth or gdatmodi.propsplt:
                gdatmodi.thislpau *= -1.
            
            gdatmodi.thislpautotl = sum(gdatmodi.thislpau)
    
    lpritotl = sum(lpri)
    
    setattr(gdatobjt, strg + 'lpritotl', lpritotl) 
    setattr(gdatobjt, strg + 'lpri', lpri)
    
    stopchro(gdat, gdatmodi, 'lpri')
    
    ### loglikelihood
    initchro(gdat, gdatmodi, 'llik')

    if not lprionly:
        
        # process a sample vector and the occupancy list to calculate secondary variables
        if gdat.elemtype == 'lens':
            lgalsour = sampvarb[getattr(gdat, strgmodl + 'indxfixplgalsour')]
            bgalsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpbgalsour')]
            specsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpspecsour')]
            sizesour = sampvarb[getattr(gdat, strgmodl + 'indxfixpsizesour')]
            ellpsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpellpsour')]
            anglsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpanglsour')]
            lgalhost = sampvarb[getattr(gdat, strgmodl + 'indxfixplgalhost')]
            bgalhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpbgalhost')]
            spechost = sampvarb[getattr(gdat, strgmodl + 'indxfixpspechost')]
            sizehost = sampvarb[getattr(gdat, strgmodl + 'indxfixpsizehost')]
            if raww:
                beinhost = 0.
            else:
                beinhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpbeinhost')]
            ellphost = sampvarb[getattr(gdat, strgmodl + 'indxfixpellphost')]
            anglhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpanglhost')]
            if raww:
                sherhost = 0.
            else:
                sherhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpsherhost')]
            sanghost = sampvarb[getattr(gdat, strgmodl + 'indxfixpsanghost')]
           
            ## host halo deflection
            defl = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost)
            
            ## external shear
            deflextr = retr_deflextr(gdat, sherhost, sanghost)
            defl += deflextr
        
            ## PSF kernel
            psfnkern = []
            for i in gdat.indxener:
                psfnkern.append(AiryDisk2DKernel(psfp[i] / gdat.sizepixl))
            
        if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
            ## PSF off-axis factor
            if oaxitype:
                onor = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfponor')]
                oind = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfpoind')]
                factoaxi = retr_factoaxi(gdat, gdat.binsoaxi, onor, oind)
        
            psfntype = getattr(gdat, strgmodl + 'psfntype')
            
            psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
            
            if oaxitype:
                psfnintp = []
                for p in gdat.indxoaxi:
                    psfnintp.append(interp1d_pick(gdat.binsanglfull, psfn[:, :, :, p], axis=1))
            else:
                psfnintp = interp1d_pick(gdat.binsanglfull, psfn, axis=1)
            setattr(gdatobjt, strg + 'psfnintp', psfnintp)
        
        if gdat.elemtype == 'lens':

            ## subhalos
            if numbpntsconc > 0 and not raww:
                if gdat.jitt:
                    timefunc('retr_deflelem_jitt', retr_deflelem_jitt, gdat.deflelem, gdat.indxpixl, gdat.lgalgrid, \
                                      gdat.bgalgrid, numbpntsconc, dicttemp['lgalconc'], dicttemp['bgalconc'], dicttemp['defsconc'], dicttemp['ascaconc'], dicttemp['acutconc'])
                    retr_deflelem_jitt(gdat.deflelem, gdat.indxpixl, gdat.lgalgrid, gdat.bgalgrid, numbpntsconc, dicttemp['lgalconc'], dicttemp['bgalconc'], \
                                                                dicttemp['defsconc'], dicttemp['ascaconc'], dicttemp['acutconc'])
                else:

                    if gdat.propwithsing and gdat.pertmodleval and strg == 'next':
                        
                        if False:
                            print 'gdat.fittnumbtrapcumr'
                            print gdat.fittnumbtrapcumr
                            print 'gdat.fittindxsampcompinit'
                            print gdat.fittindxsampcompinit
                            print 'gdatmodi.thisindxsampcomp'
                            print gdatmodi.thisindxsampcomp
                        if gdatmodi.propwith:
                            if gdatmodi.indxsampmodi < gdat.fittindxsampcompinit:
                                deflelem = gdatmodi.thisdeflelem
                                numbsubhtemp = 0
                            else:
                                deflelem = copy(gdatmodi.thisdeflelem)
                                numbsubhtemp = 2
                                for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                                    thiscomp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi]]
                                    nextcomp = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi]]
                                    if namecomp == gdat.namefeatsign:
                                        dicttemp[namecomp + 'eval'] = array([-thiscomp, nextcomp])
                                    else:
                                        dicttemp[namecomp + 'eval'] = array([thiscomp, nextcomp])
                                    
                                    if False:
                                        print 'thiscomp'
                                        print thiscomp
                                        print 'nextcomp'
                                        print nextcomp
                                        print 'dicttemp[namecomp + eval]'
                                        print dicttemp[namecomp + 'eval']
                        elif gdatmodi.propbrth or gdatmodi.propdeth:
                            deflelem = copy(gdatmodi.thisdeflelem)
                            numbsubhtemp = 1
                            if gdatmodi.propbrth:
                                comp = gdatmodi.nextsampvarb[gdatmodi.indxsamptran[0]]
                            else:
                                comp = gdatmodi.thissampvarb[gdatmodi.indxsamptran[0]]
                            for k, namecomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
                                dicttemp[namecomp + 'eval'] = empty(1)
                                if namecomp == gdat.namefeatsign and gdatmodi.propdeth:
                                    dicttemp[namecomp + 'eval'][0] = -comp[k]
                                else:
                                    dicttemp[namecomp + 'eval'][0] = comp[k]
                    else:
                        deflelem = zeros((gdat.numbpixl, 2))
                        numbsubhtemp = numbpntsconc
                        for namecomp in getattr(gdat, strgmodl + 'liststrgcomptotl'):
                            dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']
                
                    for k in range(numbsubhtemp):
                        if False and strg == 'next' and gdat.propwithsing and gdat.pertmodleval:
                            print 'k'
                            print k
                            for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                                print 'namecomp'
                                print namecomp
                                print 'comp'
                                print dicttemp[namecomp + 'eval'][k]
                        
                        ampleval = dicttemp[gdat.namecompampl + 'eval'][k]
                        indxpixltemp = retr_indxpixleval(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], ampleval, gdat.evalcirc)

                        if gdat.variasca:
                            asca = dicttemp['ascaeval'][k]
                        else:
                            asca = gdat.ascaglob

                        if gdat.variacut:
                            acut = dicttemp['acuteval'][k]
                        else:
                            acut = gdat.acutglob

                        deflelem[indxpixltemp, :] += retr_defl_jitt(gdat.indxpixl, gdat.lgalgrid, gdat.bgalgrid, \
                                                     dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], dicttemp['defseval'][k], 0., 0., \
                                                                                                                    asca=asca, acut=acut, indxpixltemp=indxpixltemp)
                    defl += deflelem
                    
                    if False and strg == 'next' and gdat.propwithsing and gdat.pertmodleval:
                        print
                        print
                        print
                        print
                        print
                        print
            
                setattr(gdatobjt, strg + 'deflelem', deflelem)

            defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
            
            # lensed image
            lensflux = retr_mapssers(gdat, gdat.lgalgridcart - defl[:, :, 0], gdat.bgalgridcart - defl[:, :, 1], lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
            
            # host emission
            hostfluxmaps = retr_mapssers(gdat, gdat.lgalgridcart, gdat.bgalgridcart, lgalhost, bgalhost, spechost, sizehost, ellphost, anglhost)
           
            # total emission
            backfluxcart = getattr(gdat, strgmodl + 'backfluxcart')
            modlfluxuncv = lensflux + bacp * backfluxcart[0] + hostfluxmaps
            
            # convolve the lensed image with the PSF
            # temp
            # modlflux = modlfluxuncv.reshape((gdat.numbener, gdat.numbpixl, 1))
            modlflux = empty_like(gdat.expo)
            for i in gdat.indxener:
                modlflux[i, :, 0] = convolve_fft(modlfluxuncv[i, :, :, 0], psfnkern[i]).flatten()
            
            if False:
                print 'modlfluxuncv'
                summgene(modlfluxuncv)
                print 'modlflux'
                summgene(modlflux)
                print
                print
                print

            setattr(gdatobjt, strg + 'defl', defl)
            
        if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
            if gdat.propwithsing and gdat.pertmodleval and strg == 'next':
                if gdatmodi.propwith:
                    if gdatmodi.indxsampmodi < gdat.fittindxsampcompinit:
                        deflelem = gdatmodi.thispntsflux
                        numbpntstemp = 0
                    else:
                        pntsflux = copy(gdatmodi.thispntsflux)
                        numbpntstemp = 2
                        for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                            thiscomp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi]]
                            nextcomp = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi]]
                            if namecomp == gdat.namefeatsign:
                                dicttemp[namecomp + 'eval'] = array([-thiscomp, nextcomp])
                            else:
                                dicttemp[namecomp + 'eval'] = array([thiscomp, nextcomp])
                            
                            if False:
                                print 'thiscomp'
                                print thiscomp
                                print 'nextcomp'
                                print nextcomp
                                print 'dicttemp[namecomp + eval]'
                                print dicttemp[namecomp + 'eval']
                elif gdatmodi.propbrth or gdatmodi.propdeth:
                    pntsflux = copy(gdatmodi.thispntsflux)
                    numbpntstemp = 1
                    if gdatmodi.propbrth:
                        comp = gdatmodi.nextsampvarb[gdatmodi.indxsamptran[0]]
                    else:
                        comp = gdatmodi.thissampvarb[gdatmodi.indxsamptran[0]]
                    for k, namecomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
                        dicttemp[namecomp + 'eval'] = empty(1)
                        if namecomp == gdat.namefeatsign and gdatmodi.propdeth:
                            dicttemp[namecomp + 'eval'][0] = -comp[k]
                        else:
                            dicttemp[namecomp + 'eval'][0] = comp[k]
            else:
                pntsflux = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
                numbpntstemp = numbpntsconc
                for namecomp in ['lgal', 'bgal', 'spec']:#getattr(gdat, strgmodl + 'liststrgcomptotl'):
                    dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']
            
            if gdat.elemtype == 'lght':
                specconctemp = dicttemp['speceval']
            else:
                specconctemp = dicttemp['nobjeval'][None, :]

            for k in range(numbpntstemp):
                indxpixleval = retr_indxpixleval(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], specconctemp, gdat.evalcirc)
                pntsflux[:, indxpixleval, :] += retr_pntsflux(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], specconctemp[:, k], psfnintp, oaxitype, indxpixleval)
        
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
       
            proc_datacnts(gdat)
        
        ### log-likelihood
        #llik = retr_llik_mult(gdat, modlcnts)
        llik = retr_llik_depr(gdat, modlcnts)
        lliktotl = sum(llik)
    
        setattr(gdatobjt, strg + 'llik', llik) 
    else:
        lliktotl = 0.

    setattr(gdatobjt, strg + 'lliktotl', lliktotl) 
    
    stopchro(gdat, gdatmodi, 'llik')

    if lprionly:
        return
    else:
        lpostotl = lpritotl + lliktotl
        setattr(gdatobjt, strg + 'lpostotl', lpostotl) 

    if strg == 'next':
        setattr(gdatmodi, 'thislpriprop', lpri)
    
    stopchro(gdat, gdatmodi, 'proc')
    
    if fast:
        return

    ## tertiary variables that are not needed for likelihood evaluation
    if strg != 'next':
        
        liststrgfeatcorr = getattr(gdat, strgmodl + 'liststrgfeatcorr')
        liststrgfeatcorrtotl = getattr(gdat, strgmodl + 'liststrgfeatcorrtotl')
       
        setattr(gdatobjt, strg + 'bacp', bacp)
        setattr(gdatobjt, strg + 'psfp', psfp)
        resicnts = gdat.datacnts - modlcnts
        setattr(gdatobjt, strg + 'resicnts', resicnts)
        
        ## derived variables
        if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
            if oaxitype:
                setattr(gdatobjt, strg + 'factoaxi', factoaxi)
           
            setattr(gdatobjt, strg + 'psfn', psfn)

            ### PSF FWHM
            fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            setattr(gdatobjt, strg + 'fwhm', fwhm)
            
	    	### mean PS flux map 
            pntsfluxmean = sum(sum(pntsflux * gdat.expo, 2), 1) / sum(sum(gdat.expo, 2), 1)
            setattr(gdatobjt, strg + 'pntsfluxmean', pntsfluxmean)

            if gdat.calcerrr and gdat.fittnumbtrap > 0:
                if False:
                    pntsfluxfull = retr_pntsfluxtotl(gdat, dicttemp['lgalconc'], dicttemp['bgalconc'], dicttemp['specconc'], psfnintp, gdat.fittoaxitype, evalcirc=False)
                else:
                    pntsfluxfull = zeros_like(gdat.datacnts)
                pntscnts = retr_cntsmaps(gdat, pntsflux)
                pntscntsfull = retr_cntsmaps(gdat, pntsfluxfull)
                errrcnts = pntscnts - pntscntsfull
                setattr(gdatobjt, strg + 'errrcnts', errrcnts)
                if False and amax(fabs(errrcnts)) > 0.1:
                    raise Exception('Approximation error in calculating the PS flux map is above the tolerance level.')

            #fluxbrgt, fluxbrgtassc = retr_fluxbrgt(gdat, dicttemp['lgalconc'], dicttemp['bgalconc'], concatenate(dicttemp['flux']))
            #setattr(gdatobjt, strg + 'fluxbrgt', fluxbrgt)
            #setattr(gdatobjt, strg + 'fluxbrgtassc', fluxbrgtassc)

        if gdat.elemtype == 'lens':
            modlcntsuncv = retr_cntsmaps(gdat, modlfluxuncv, cart=True)
            lenscnts = retr_cntsmaps(gdat, lensflux, cart=True)
            hostcntsmaps = retr_cntsmaps(gdat, hostfluxmaps, cart=True)
            
            setattr(gdatobjt, strg + 'psfnkern', psfnkern)
            setattr(gdatobjt, strg + 'modlcntsuncv', modlcntsuncv)
            setattr(gdatobjt, strg + 'lenscnts', lenscnts)
            setattr(gdatobjt, strg + 'hostcntsmaps', hostcntsmaps)

            ### sort with respect to deflection at scale radius
            if numbpntsconc > 0:
                indxpntssortbrgt = argsort(dicttemp[gdat.namefeatsort + 'conc'])[::-1]
                for strgcomp in liststrgcomptotl:
                    dicttemp[strgcomp + 'sort'] = dicttemp[strgcomp + 'conc'][indxpntssortbrgt][:numbpntsconc]

            ### mass budget
            if gdat.variasca:
                asca = dicttemp['asca'][0]
            else:
                asca = gdat.ascaglob
            if gdat.variacut:
                acut = dicttemp['acut'][0]
            else:
                acut = gdat.acutglob
            factmcutfromdefs = retr_factmcutfromdefs(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour, asca, acut) 
            masssubh = factmcutfromdefs * dicttemp['defs'][0]
            masssubhtotl = array([sum(masssubh)])
            masshostbein = array([gdat.massfrombein * beinhost**2])
            fracsubh = masssubhtotl / masshostbein
            
            setattr(gdatobjt, strg + 'masssubhtotl', masssubhtotl)
            setattr(gdatobjt, strg + 'masshostbein', masshostbein)
            setattr(gdatobjt, strg + 'fracsubh', fracsubh)

            deflsing = zeros((gdat.numbpixl, 2, gdat.numbdeflsingplot))
            numbdeflsing = min(gdat.numbdeflpntsplot, numbpntsconc) + 2
            if numbpntsconc > 0:
                numbdeflsing += 1
            for k in range(numbdeflsing):
                
                if False:
                    print 'k'
                    print k
                
                if k == 0:
                    deflhost = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost)
                    deflsing[:, :, k] = deflhost
                elif k == 1:
                    deflsing[:, :, k] = deflextr
                elif k == 2:
                    deflsing[:, :, k] = defl.reshape((gdat.numbpixl, 2)) - deflextr - deflhost
                else:
                    if gdat.evalcirc == 'full':
                        indxpixltemp = gdat.indxpixl
                    else:
                        indxpixlpnts = retr_indxpixl(gdat, dicttemp['bgalsort'][k-2], dicttemp['lgalsort'][k-2])
                        indxproxtemp = digitize(dicttemp['defssort'][k-2], gdat.binsprox) - 1
                        indxpixltemp = gdat.indxpixlprox[indxproxtemp][indxpixlpnts]
                        if isinstance(indxpixltemp, int):
                            indxpixltemp = gdat.indxpixl
                    
                    if gdat.variasca:
                        asca = dicttemp['ascasort'][k-3]
                    else:
                        asca = gdat.ascaglob

                    if gdat.variacut:
                        acut = dicttemp['acutsort'][k-3]
                    else:
                        acut = gdat.acutglob

                    deflsing[indxpixltemp, :, k] = retr_defl(gdat, dicttemp['lgalsort'][k-3], dicttemp['bgalsort'][k-3], dicttemp['defssort'][k-3], 0., 0., \
                                                                                                                    asca=asca, acut=acut, indxpixltemp=indxpixltemp)
                    
            deflsing = deflsing.reshape((gdat.numbsidecart, gdat.numbsidecart, 2, gdat.numbdeflsingplot))

            ### convergence
            deflelem = deflelem.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
        
            magn = 1. / retr_invm(gdat, defl) 
            conv = retr_conv(gdat, defl) 
            convelem = retr_conv(gdat, deflelem) 
            
            convpsec = retr_psec(gdat, conv)
            convpsecelem = retr_psec(gdat, convelem)
            convpsecodim = retr_psecodim(gdat, convpsec) 
            convpsecelemodim = retr_psecodim(gdat, convpsecelem) 
            
            print 'defl'
            print defl
            print 'gdat.binsdefl'
            print gdat.binsdefl
            print

            histdefl = histogram(defl, bins=gdat.binsdefl)[0]
            histdeflelem = histogram(deflelem, bins=gdat.binsdeflelem)[0]
            
            setattr(gdatobjt, strg + 'magn', magn)
            setattr(gdatobjt, strg + 'conv', conv)
            setattr(gdatobjt, strg + 'convelem', convelem)
            setattr(gdatobjt, strg + 'convpsec', convpsec)
            setattr(gdatobjt, strg + 'convpsecelem', convpsecelem)
            setattr(gdatobjt, strg + 'convpsecodim', convpsecodim)
            setattr(gdatobjt, strg + 'convpsecelemodim', convpsecelemodim)
            
            setattr(gdatobjt, strg + 'histdefl', histdefl)
            setattr(gdatobjt, strg + 'histdeflelem', histdeflelem)
            setattr(gdatobjt, strg + 'deflsing', deflsing)
        
        liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
    
        if gdat.elemtype == 'lens':
            lensfluxmean = mean(sum(lensflux, 3), 0)
            
        ## element features
        if numbtrap > 0:
            
            ### derived quantities
            for l in range(numbpopl):
                #### radial and angular coordinates
                dicttemp['gang'][l] = retr_gang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                dicttemp['aang'][l] = retr_aang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                
                if gdat.elemtype == 'lght':
                    #### number of expected counts
                    dicttemp['cnts'][l] = retr_pntscnts(gdat, dicttemp['lgal'][l], dicttemp['bgal'][l], dicttemp['spec'][l])
                
            #### delta log-likelihood
            gdatmoditemp = tdpy.util.gdatstrt()
            for l in range(numbpopl):
                dicttemp['deltllik'][l] = zeros(numbpnts[l])
                for k in range(numbpnts[l]):
                    if gdat.elemtype == 'lens':
                        if gdat.variasca:
                            asca = dicttemp['asca'][l][k]
                        else:
                            asca = gdat.ascaglob

                        if gdat.variacut:
                            acut = dicttemp['acut'][l][k]
                        else:
                            acut = gdat.acutglob
                        
                        defltemp = copy(defl.reshape((gdat.numbpixl, 2)))

                        defltemp[indxpixltemp, :] -= retr_defl(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['defs'][l][k], \
                                                                                 0., 0., asca=asca, acut=acut, indxpixltemp=indxpixltemp)
                        defltemp = defltemp.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
                        lensflux = retr_mapssers(gdat, gdat.lgalgridcart - defltemp[:, :, 0], gdat.bgalgridcart - defltemp[:, :, 1], \
                                                                                                            lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)
                        modlfluxuncv = lensflux + bacp * backfluxcart[0] + hostfluxmaps
                        modlflux = empty_like(gdat.expo)
                        for i in gdat.indxener:
                            modlflux[i, :, 0] = convolve_fft(modlfluxuncv[i, :, :, 0], psfnkern[i]).flatten()
                    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
                        pntsfluxtemp = copy(pntsflux)

                        if gdat.elemtype == 'lght':
                            spectemp = dicttemp['spec'][l][:, k]
                        else:
                            spectemp = dicttemp['nobj'][l][None, k]
                        pntsfluxtemp[:, indxpixleval, :] -= retr_pntsflux(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], spectemp, psfnintp, oaxitype, indxpixleval)
                        modlflux = retr_mapslght(gdat, bacp, pntsfluxtemp, gdat.indxcube)
                    modlcnts = retr_cntsmaps(gdat, modlflux)
                    #nextllik = retr_llik_mult(gdat, modlcnts)
                    nextllik = retr_llik_depr(gdat, modlcnts)
                    nextlliktotl = sum(nextllik)
                    dicttemp['deltllik'][l][k] = lliktotl - nextlliktotl
            
            if gdat.elemtype == 'lght':
                #### spectra
                for l in indxpopl:
                    dicttemp['specplot'][l] = retr_spec(gdat, dicttemp['flux'][l], dicttemp['sind'][l], dicttemp['curv'][l], dicttemp['expo'][l], spectype[l], plot=True)
                
            if gdat.elemtype == 'lens':
                #### distance to the source
                for l in range(numbpopl):
                    dicttemp['diss'][l] = retr_angldist(gdat, dicttemp['lgal'][l],  dicttemp['bgal'][l], lgalsour, bgalsour)
                
                for l in indxpopl:
                    dicttemp['deflprof'][l] = empty((gdat.numbangl, numbpnts[l]))
                    dicttemp['mcut'][l] = empty(numbpnts[l])
                    dicttemp['dotn'][l] = empty(numbpnts[l])
                    dicttemp['dots'][l] = empty(numbpnts[l])
                    dicttemp['dotv'][l] = empty(numbpnts[l])
                    dicttemp['dotm'][l] = empty(numbpnts[l])
                    deflsubh = zeros((gdat.numbpixl, 2, numbpnts[l]))
                    for k in arange(numbpnts[l]):
                        
                        if gdat.variasca:
                            asca = dicttemp['asca'][l][k]
                        else:
                            asca = gdat.ascaglob

                        if gdat.variacut:
                            acut = dicttemp['acut'][l][k]
                        else:
                            acut = gdat.acutglob
                        
                        #### deflection profiles
                        dicttemp['deflprof'][l][:, k] = retr_deflcutf(gdat.binsangl[1:], dicttemp['defs'][l][k], asca, acut)
             
                        ### truncated mass 
                        dicttemp['mcut'][l][k] = retr_mcut(gdat, dicttemp['defs'][l][k], asca, acut)
                
                        #### dot product with the source flux gradient
                        # temp -- weigh the energy and PSF bins
                        dicttemp['dotn'][l][k] = retr_dots(gdat, gdat.datacntscart[0, :, :, 0], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], 1., \
                                                                                                 asca, acut, gdat.indxpixl)
                        dicttemp['dots'][l][k] = dicttemp['defs'][l][k] * dicttemp['dotn'][l][k] / gdat.sizepixl
                        dicttemp['dotv'][l][k] = retr_dots(gdat, gdat.datacntscart[0, :, :, 0], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], 1., \
                                                                                                 asca, acut, gdat.indxpixl, absv=True)
                        dicttemp['dotm'][l][k] = retr_dots(gdat, lenscnts[0, :, :, 0], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], 1., \
                                                                                                 asca, acut, gdat.indxpixl)
                   
            if strgmodl == 'true':
                indx = [[] for l in indxpopl]
                for l in indxpopl:
                    indx[l] = where(dicttemp['deltllik'][l] > 0.5 * numbcomp[l])[0]
                    dicttemp[gdat.namefeatsign + 'sign'][l] = zeros(numbpnts[l])
                    dicttemp[gdat.namefeatsign + 'sign'][l][indx[l]] = dicttemp[gdat.namefeatsign][l][indx[l]]
                setattr(gdat, 'trueindxelemsign', indx)
                
            # temp
            if False and strgmodl == 'true':
                for strgfeat in liststrgfeattotl:
                    minm = getattr(gdat, 'minm' + strgfeat)
                    maxm = getattr(gdat, 'maxm' + strgfeat)
                    if strgfeat == 'gang':
                        minm *= gdat.anglfact
                        maxm *= gdat.anglfact
                        dicttemp[strgfeat][l] *= gdat.anglfact
                    for l in indxpopl:
                        if where(minm > dicttemp[strgfeat][l])[0].size > 0 or where(maxm < dicttemp[strgfeat][l])[0].size > 0:
                            print 'Element feature outside the limits'
                            print 'strgfeat'
                            print strgfeat
                            print 'minm'
                            print minm
                            print 'maxm'
                            print maxm
                            print 'dicttemp[strgfeat][l]'
                            print dicttemp[strgfeat][l]
                            print

            ### distribution of element parameters and features
            #### find the indices of the model PSs that are in the comparison area
            # temp -- this should be done for true sources
            indxmodlpntscomp = [[] for l in range(numbpopl)]
            for l in range(numbpopl):
                indxmodlpntscomp[l] = where((fabs(dicttemp['lgal'][l]) < gdat.maxmgangdata) & (fabs(dicttemp['bgal'][l]) < gdat.maxmgangdata))[0]
            
            #### one dimensional
            for strgfeat in liststrgfeattotl:
                if strgfeat == 'spec':
                    temp = zeros((numbpopl, gdat.numbbinsplot, gdat.numbener))
                elif strgfeat == 'cnts':
                    temp = zeros((numbpopl, gdat.numbbinsplot, gdat.numbener, gdat.numbevtt))
                else:
                    temp = zeros((numbpopl, gdat.numbbinsplot))
                dicttemp['hist' + strgfeat] = temp
            
                for l in indxpopl:
                    if strgfeat == 'specplot' or strgfeat == 'deflprof':
                        continue
                    elif strgfeat == 'spec':
                        for i in gdat.indxener:
                            dicttemp['hist' + strgfeat][l, :, i] = histogram(dicttemp['spec'][l][i, indxmodlpntscomp[l]], gdat.binsspec)[0]
                    elif strgfeat == 'cnts':
                        for i in gdat.indxener:
                            for m in gdat.indxevtt:
                                dicttemp['hist' + strgfeat][l, :, i, m] = histogram(dicttemp['cnts'][l][i, indxmodlpntscomp[l], m], gdat.binscnts)[0]
                    elif not (strgfeat == 'curv' and spectype[l] != 'curv' or strgfeat == 'expo' and spectype[l] != 'expo'):
                        bins = getattr(gdat, 'bins' + strgfeat)
                        if strgfeat == gdat.namefeatsign + 'sign':
                            indx = intersect1d(getattr(gdat, strgmodl + 'indxelemsign'), indxmodlpntscomp[l])
                            dicttemp['hist' + strgfeat][l, :] = histogram(dicttemp[strgfeat[:-4]][l][indx], bins)[0]
                        else:
                            dicttemp['hist' + strgfeat][l, :] = histogram(dicttemp[strgfeat][l][indxmodlpntscomp[l]], bins)[0]
                
            #### two dimensional
            for strgfrst in liststrgfeatcorrtotl:
                for strgseco in liststrgfeatcorrtotl:
                    dicttemp['hist' + strgfrst + strgseco] = zeros((numbpopl, gdat.numbbinsplot, gdat.numbbinsplot))
            
            for l in indxpopl:
                for a, strgfrst in enumerate(liststrgfeatcorr[l]):
                    for b, strgseco in enumerate(liststrgfeatcorr[l]):
                        if a < b:
                            binsfrst = getattr(gdat, 'bins' + strgfrst)
                            binsseco = getattr(gdat, 'bins' + strgseco)
                            dicttemp['hist' + strgfrst + strgseco][l, :, :] = histogram2d(dicttemp[strgfrst][l][indxmodlpntscomp[l]], \
                                                                                                        dicttemp[strgseco][l][indxmodlpntscomp[l]], [binsfrst, binsseco])[0]
                            setattr(gdatobjt, strg + 'hist' + strgfrst + strgseco, dicttemp['hist' + strgfrst + strgseco])
            
            ### priors on element parameters and features
            for strgfeat in liststrgfeattotl:
                dicttemp['hist' + strgfeat + 'prio'] = empty((numbpopl, gdat.numbbinsplotprio))
            
                for l in range(numbpopl):
                    if strgfeat in liststrgfeatprio[l]:
                        minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                        maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                        bins = getattr(gdat, 'bins' + strgfeat)
                        delt = getattr(gdat, 'delt' + strgfeat)
                        
                        xdat = getattr(gdat, strgmodl + 'mean' + strgfeat + 'prio')
                        deltprio = getattr(gdat, strgmodl + 'delt' + strgfeat + 'prio')
                        
                        booltemp = False
                        if strgfeat == 'gang' and spatdisttype[l] == 'gang' or strgfeat == 'bgal' and spatdisttype[l] == 'disc': 
                            scal = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distscal')[l]]
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
                        # temp 
                        if strgfeat == 'flux' or strgfeat == 'defs' or strgfeat == 'nobj':
                            slop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
                            if ampldisttype == 'powr':
                                pdfn = pdfn_powr(xdat, minm, maxm, slop)
                            if ampldisttype == 'igam':
                                cutf = getattr(gdat, 'cutf' + strgfeat)
                                pdfn = pdfn_igam(xdat, slop, cutf)
                            booltemp = True
                        elif strgfeat == 'sind' or strgfeat == 'curv' and spectype[l] == 'curv' or strgfeat == 'expo' and spectype[l] == 'expo' or strgfeat == 'asca' or \
                                                                                                                                strgfeat == 'acut':
                            # this does not work for mismodeling
                            meanvarb = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distmean')[l]]
                            stdv = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'diststdv')[l]]
                            if strgfeat == 'expo' and spectype[l] == 'expo':
                                pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                            else:
                                pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                            booltemp = True
                        
                        if booltemp:
                            dicttemp['hist' + strgfeat + 'prio'][l, :] = meanpnts[l] * pdfn * deltprio * delt[0] / deltprio[0]
        
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
        ## copy element features to the global object
        for strgfeat in liststrgfeattotl:
            feat = [[] for l in indxpopl]
            for l in indxpopl:
                if strgfeat in liststrgfeat[l]:
                    if strg == 'true' and not isinstance(dicttemp[strgfeat][l], list):
                        shap = list(ones(dicttemp[strgfeat][l].ndim, dtype=int))
                        feat[l] = tile(dicttemp[strgfeat][l], [3] + shap)
                    if strg == 'this':
                        feat[l] = dicttemp[strgfeat][l]
            setattr(gdatobjt, strg + strgfeat, feat)
                
        for strgfeat in liststrgfeattotl:
            setattr(gdatobjt, strg + 'hist' + strgfeat, dicttemp['hist' + strgfeat])
            setattr(gdatobjt, strg + 'hist' + strgfeat + 'prio', dicttemp['hist' + strgfeat + 'prio'])

        ### PS indices to compare with the reference catalog
        if strg == 'this' and numbtrap > 0 and gdat.trueinfo:
            
            if gdat.elemtype == 'lens' and gdat.datatype == 'mock':
                gdatmodi.thisdeflsingresi = gdatmodi.thisdeflsing - gdat.truedeflsing
                gdatmodi.thisdeflresi = gdatmodi.thisdefl - gdat.truedefl
                gdatmodi.thisconvelemresi = gdatmodi.thisconvelem - gdat.trueconvelem
                gdatmodi.thismagnresi = gdatmodi.thismagn - gdat.truemagn
                gdatmodi.thisconvelempercresi = 100. * fabs(gdatmodi.thisconvelemresi / gdat.trueconvelem)
                gdatmodi.thismagnpercresi = 100. * fabs(gdatmodi.thismagnresi / gdat.truemagn)
                
                if False:
                    print 'proc_samp'
                    print 'gdat.truelgal'
                    print gdat.truelgal
                    print 'lgal'
                    print dicttemp['lgal']
                    print 'gdat.truedefs'
                    print gdat.truedefs
                    print 'defs'
                    print dicttemp['defs']
                    print 'gdat.truedeflsing'
                    summgene(gdat.truedeflsing)
                    print 'gdat.truedefl'
                    summgene(gdat.truedefl)
                    print 'gdatmodi.thisdefl'
                    summgene(gdatmodi.thisdefl)
                    print 'gdatmodi.thisdeflresi'
                    summgene(gdatmodi.thisdeflresi)
                    for k in range(gdat.numbdeflsingplot):
                        print 'k'
                        print k
                        print 'gdatmodi.thisdeflsing'
                        summgene(gdatmodi.thisdeflsing[:, :, :, k])
                        print 'gdatmodi.thisdeflsingresi'
                        summgene(gdatmodi.thisdeflsingresi[:, :, :, k])
                    print
                    print
                    print
                
                gdatmodi.thisdeflcomp = 180. / pi * sqrt((1. - sum(gdatmodi.thisdefl * gdat.truedefl, axis=2) / sqrt(sum(gdatmodi.thisdefl**2, axis=2)) / \
                                            sqrt(sum(gdat.truedefl**2, axis=2))) / 2.)
    
    # correlate the catalog sample with the reference catalog
    # temp
    if strg == 'this' and gdat.trueinfo:
        gdatmodi.trueindxpntsasscmiss = [[] for l in gdat.trueindxpopl]
        gdatmodi.trueindxpntsasschits = [[] for l in gdat.trueindxpopl]
        featassc = dict()
        for strgfeat in liststrgfeatodimtotl:
            featassc[strgfeat] = [[] for l in gdat.trueindxpopl]
        cmpl = [[] for l in gdat.trueindxpopl]
        indxfittpntsassc = [[] for l in gdat.trueindxpopl]
        indxfittpntsfals = [[] for l in gdat.trueindxpopl]
        fdis = [[] for l in gdat.trueindxpopl]

        for l in gdat.trueindxpopl:
            indxmodlpnts = zeros(gdat.truenumbpnts[l], dtype=int) - 1
            numbassc = zeros(gdat.truenumbpnts[l])
            metrassc = zeros(gdat.truenumbpnts[l]) + 3 * gdat.maxmgang
            
            for k in range(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[l]].astype(int)):
               
                # determine which true elements satisfy the match criterion
                metr = retr_angldist(gdat, gdat.truelgal[l][0, :], gdat.truebgal[l][0, :], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k])

                trueindxpntstemp = where(metr < gdat.anglassc)[0]
                
                if trueindxpntstemp.size > 0:
                    # if there are multiple associated true elements, sort them
                    indx = argsort(metr[trueindxpntstemp])
                    metr = metr[trueindxpntstemp][indx]
                    trueindxpntstemp = trueindxpntstemp[indx]
                    
                    # store the index of the model PS
                    numbassc[trueindxpntstemp[0]] += 1
                    if metr[0] < metrassc[trueindxpntstemp[0]]:
                        metrassc[trueindxpntstemp[0]] = metr[0]
                        indxmodlpnts[trueindxpntstemp[0]] = k
        
            # collect the associated element features
            indx = where(indxmodlpnts >= 0)[0]
            for strgfeat in liststrgfeatodim[l]:
                featassc[strgfeat][l] = zeros(gdat.truenumbpnts[l])
                featassc[strgfeat][l][indx] = dicttemp[strgfeat][l][indxmodlpnts[indx]]
                
            # divide associations into subgroups
            for k in range(gdat.truenumbpnts[l]):
                if numbassc[k] == 0:
                    gdatmodi.trueindxpntsasscmiss[l].append(k)
                else:
                    gdatmodi.trueindxpntsasschits[l].append(k)
            temp = where(indxmodlpnts >= 0)[0]
            featassc[strgfeat][l][temp] = dicttemp[strgfeat][l][indxmodlpnts[temp]]

            indxfittpntsassc[l] = unique(indxmodlpnts[where(indxmodlpnts > -1)])
            indxfittpntsfals[l] = setdiff1d(arange(numbpnts[l]), indxfittpntsassc[l])
            
            cmpl[l] = array([float(len(gdatmodi.trueindxpntsasschits[l])) / gdat.truenumbpnts[l]])
            fdis[l] = array([float(indxfittpntsfals[l].size) / numbpnts[l]])
            
            if False:
                print 'gdatmodi.trueindxpntsasschits[l]'
                print gdatmodi.trueindxpntsasschits[l]
                print 'gdat.truenumbpnts[l]'
                print gdat.truenumbpnts[l]
                print 'cmpl'
                print cmpl
                print 'indxfittpntsfals[l]'
                print indxfittpntsfals[l]
                print 'numbpnts[l]'
                print numbpnts[l]
                print 'fdis'
                print fdis
                print

        for l in gdat.trueindxpopl:
            setattr(gdatmodi, 'thiscmplpop%d' % l, cmpl[l])
        
        for l in gdat.fittindxpopl:
            setattr(gdatmodi, 'thisfdispop%d' % l, fdis[l])
        
        for strgfeat in liststrgfeatodimtotl:
            setattr(gdatmodi, 'this' + strgfeat + 'assc', featassc[strgfeat])
       
        # completeness
        for strgfeat in liststrgfeatodimtotl:
            cmplfeat = empty((gdat.truenumbpopl, gdat.numbbinsplot))
            truehistfeat = getattr(gdat, 'truehist' + strgfeat)
            for l in gdat.trueindxpopl:
                indx = where(isfinite(featassc[strgfeat][l]))[0]
                bins = getattr(gdat, 'bins' + strgfeat)
                print 'strgfeat'
                print strgfeat
                print 'bins'
                print bins
                print 'dicttemp[strgfeat]'
                print dicttemp[strgfeat]
                print 'indxfittpntsassc'
                print indxfittpntsassc
                print

                histfeatassc = histogram(dicttemp[strgfeat][l][indxfittpntsassc[l]], bins=bins)[0]
                cmplfeat[l, :] = histfeatassc / truehistfeat[l, :]
            setattr(gdatmodi, 'thiscmpl' + strgfeat, cmplfeat)
            
        # false discovery rate
        for strgfeat in liststrgfeatodimtotl:
            fdisfeat = empty((gdat.fittnumbpopl, gdat.numbbinsplot))
            for l in gdat.fittindxpopl:
                indx = where(isfinite(featassc[strgfeat][l]))[0]
                bins = getattr(gdat, 'bins' + strgfeat)
                histfeatfals = histogram(dicttemp[strgfeat][l][indxfittpntsfals[l]], bins=bins)[0]
                fitthistfeat = getattr(gdatmodi, 'thishist' + strgfeat)
                fdisfeat[l, :] = histfeatfals / fitthistfeat[l, :]
            setattr(gdatmodi, 'thisfdis' + strgfeat, fdisfeat)


def proc_datacnts(gdat):

    if gdat.pixltype == 'cart':
        gdat.datacntscart = gdat.datacnts.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
    gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix
    if gdat.enerdiff:
        gdat.datafluxmean /= gdat.deltener


def retr_info(pdfnpost, pdfnprio):
    
    info = pdfnpost * log(pdfnpost / pdfnprio)

    return info


def retr_llik_depr(gdat, modlcnts):
    
    if gdat.liketype == 'pois':
    	llik = gdat.datacnts * log(modlcnts) - modlcnts
    if gdat.liketype == 'gaus':
        llik = -0.5 * (gdat.datacnts - modlcnts)**2 / gdat.datacnts
    
    return llik


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
   
    sersprof = spec / pi / size**2 * exp(-(1.992 * indx - 0.3271) * ((angl / size)**(1. / indx) - 1.))
     
    return sersprof



