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


def retr_indxsampcomp(gdat, indxelemfull, strgmodl):
    
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
        indxsamptemp = indxsampcompinit + numbtrapcuml[l] + array(indxelemfull[l], dtype=int) * numbcomp[l]
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
        indxsampcomp['comp'].append(repeat(indxsamptemp, numbcomp[l]) + tile(indxcomp[l], len(indxelemfull[l])))
             
    return indxsampcomp


def retr_plotpath(gdat, gdatmodi, strg, strgplot, nameinte=''):
    
    if strg == 'true' or strg == '':
        path = gdat.pathinit + nameinte + strgplot + '.pdf'
    elif strg == 'post' or strg == 'mlik':
        path = gdat.pathplot + gdat.namesampdist + '/finl/' + nameinte + strg + strgplot + '.pdf'
    elif strg == 'this':
        path = gdat.pathplot + gdat.namesampdist + '/fram/' + nameinte + strg + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


def retr_pntsfluxtotl(gdat, lgal, bgal, spec, psfnintp, oaxitype, evalcirc):

    pntsflux = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    numbpnts = lgal.size
    for k in range(numbpnts):
        indxpixleval = retr_indxpixleval(gdat, lgal[k], bgal[k], spec[gdat.indxenerfluxdist[0], k], evalcirc)
        pntsflux[:, indxpixleval, :] += retr_pntsflux(gdat, lgal, bgal, spec, psfnintp, oaxitype, indxpixleval)
    

def retr_indxpixleval(gdat, lgal, bgal, ampl, evalcirc):
    
    if not isscalar(lgal):
        lgal = lgal[0]
        bgal = bgal[0]
        ampl = ampl[0]
     
    if evalcirc == 'locl':  
        indxpixlpnts = retr_indxpixl(gdat, bgal, lgal)
        indxfluxproxtemp = digitize(ampl, gdat.binsprox) - 1

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
    
    if not isfinite(icdf).all():
        print 'icdf_igam'
        print 'icdf'
        print icdf
        print 'slop'
        print slop
        print 'cutf'
        print cutf
        print 'xdat'
        print xdat
        print
    
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
    gdatmodi.propfixp = False
    gdatmodi.propwith = False
    gdatmodi.propbrth = False
    gdatmodi.propdeth = False
    gdatmodi.propsplt = False
    gdatmodi.propmerg = False
    
    gdatmodi.thisljcbfact = 0.
    gdatmodi.thislpautotl = 0. 
    gdatmodi.thislcomfact = 0.
   
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
                gdatmodi.propsplt = True
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
            if gdat.propfixp and gdat.propcomp:
                gdatmodi.thisindxsampfull = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp['comp'])))
            elif gdat.propcomp and not gdat.propfixp:
                gdatmodi.thisindxsampfull = concatenate(gdatmodi.thisindxsampcomp['comp'])
            elif not gdat.propcomp and gdat.propfixp:
                gdatmodi.thisindxsampfull = gdat.indxfixpprop
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
    
    if gdat.diagmode:
        if gdat.probbrde == 0. and (gdatmodi.propbrth or gdatmodi.propdeth):
            raise Exception('')


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
    
        
def retr_sampvarb(gdat, strgmodl, samp, indxsampcomp=None):
    
    sampvarb = zeros_like(samp)
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    
    for k in indxfixp:
        sampvarb[k] = icdf_fixp(gdat, strgmodl, samp[k], k)
    
    if indxsampcomp != None:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
        listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
        retr_sampvarbcomp(gdat, strgmodl, indxsampcomp, indxpopl, liststrgcomp, listscalcomp, samp, sampvarb)

    return sampvarb
    

def retr_sampvarbcomp(gdat, strgmodl, indxsampcomp, indxpopl, liststrgcomp, listscalcomp, samp, sampvarb, indxelemfull=slice(None)):
    
    for l in indxpopl:
        for k, strgcomp in enumerate(liststrgcomp[l]):
            
            if listscalcomp[l][k] == 'self' or listscalcomp[l][k] == 'dexp' or listscalcomp[l][k] == 'expo' or listscalcomp[l][k] == 'powrslop':
                minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
                if listscalcomp[l][k] == 'powrslop':
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


def cdfn_comp(gdat, gdatmodi, comp, init=True):
    
    listscalcomp = gdat.fittlistscalcomp[gdatmodi.indxpoplmodi]
    cdfn = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
    for k, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
        
        if init and k > 2:
            break
            
        if listscalcomp[k] == 'self' or listscalcomp[k] == 'dexp' or listscalcomp[k] == 'expo' or listscalcomp[k] == 'powrslop':
            minm = getattr(gdat, 'fittminm' + strgcomp)
            if listscalcomp[k] == 'powrslop':
                maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                distslop = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[gdatmodi.indxpoplmodi]]
                cdfn[k] = cdfn_powr(comp[k], minm, maxm, distslop)
            else:
                fact = getattr(gdat, 'fittfact' + strgcomp)
                cdfn[k] = cdfn_self(comp[k], minm, fact)
        if listscalcomp[k] == 'igam':
            distslop = gdatmodi.sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[gdatmodi.indxpoplmodi]]
            cutf = getattr(gdat, 'cutf' + strgcomp)
            cdfn[k] = cdfn_igam(comp[k], distslop, cutf)
        if listscalcomp[k] == 'gaus':
            distmean = gdatmodi.sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[gdatmodi.indxpoplmodi]]
            diststdv = gdatmodi.sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[gdatmodi.indxpoplmodi]]
            cdfn[k] = cdfn_gaus(comp[k], distmean, diststdv)
    
    return cdfn


def retr_mrkrsize(gdat, compampl):
    
    minm = getattr(gdat, 'minm' + gdat.namecompampl) 
    maxm = getattr(gdat, 'maxm' + gdat.namecompampl) 
    mrkrsize = (sqrt(compampl) - sqrt(minm)) / (sqrt(maxm) - sqrt(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
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
    
    parastrg = ['score', 'gcore', 'stail', 'gtail', 'ntail']
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
            fermform[i, m, 4] = 1. / (1. + fermform[i, m, 4] * fermform[i, m, 2]**2 / fermform[i, m, 0]**2)

    # calculate the scale factor
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # store the fermi PSF parameters
    gdat.exprpsfp = zeros(gdat.numbener * gdat.numbevtt * numbpsfpform)
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            #if k == 0 or k == 2:
            #    gdat.exprpsfp[indxfermpsfptemp] = fermform[:, m, k] * gdat.fermscalfact[:, m]
            #else:
            #    gdat.exprpsfp[indxfermpsfptemp] = fermform[:, m, k]
            gdat.exprpsfp[indxfermpsfptemp] = fermform[:, m, k]
    

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
    
    if gdat.fittnumbtrap > 0:
        gdatmodi.thisindxelemfull = deepcopy(gdatmodi.nextindxelemfull)
        gdatmodi.thisindxsampcomp = deepcopy(gdatmodi.nextindxsampcomp)
    
        if gdat.calcllik and gdat.propwithsing:
            if gdat.elemtype == 'lght':
                gdatmodi.thispntsflux = copy(gdatmodi.nextpntsflux)
            if gdat.elemtype == 'lens':
                gdatmodi.thisdeflelem = copy(gdatmodi.nextdeflelem)


def retr_listpair(gdat, lgal, bgal):
    
    if gdat.verbtype > 1:
        print 'Finding element pairs inside the linking length...'
    
    listpair = []
    for k in range(lgal.size):
        # temp -- linking uses the Cartesian approximation, which is accurate enough for splits and merges inside a small circle
        #indxpnts = k + 1 + where(sqrt((bgal[k+1:] - bgal[k])**2 + (lgal[k+1:] - lgal[k])**2) < gdat.radispmr)[0]
        indxpnts = k + 1 + where((fabs(bgal[k+1:] - bgal[k]) < gdat.radispmr) & (fabs(lgal[k+1:] - lgal[k]) < gdat.radispmr))[0]
        for n in range(indxpnts.size):
            listpair.append([k, indxpnts[n]])
    
    if gdat.verbtype > 1:
        numbpair = len(listpair)
        print '%d pairs found.' % numbpair

    if gdat.diagmode:
        boolgood = True
        for n in range(len(listpair)):
            dist = sqrt((lgal[listpair[n][0]] - lgal[listpair[n][1]])**2 + (bgal[listpair[n][0]] - bgal[listpair[n][1]])**2)
            
            if gdat.verbtype > 1:
                print 'n'
                print n
                print 'lgal[listpair[n][0]]'
                print lgal[listpair[n][0]] * gdat.anglfact
                print 'lgal[listpair[n][1]]'
                print lgal[listpair[n][1]] * gdat.anglfact
                print 'bgal[listpair[n][0]]'
                print bgal[listpair[n][0]] * gdat.anglfact
                print 'bgal[listpair[n][1]]'
                print bgal[listpair[n][1]] * gdat.anglfact
                print 'dist'
                print dist * gdat.anglfact

            if dist >= gdat.radispmr:
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
    

    gdat.exprspec = zeros((3, 2, gdat.exprlgal.size))
    gdat.exprspec[0, 0, :] = fluxchansoft * 0.624e9
    gdat.exprspec[0, 1, :] = fluxchanhard * 0.624e9 / 16.
    # temp
    gdat.exprspec[1, :, :] = gdat.exprspec[0, :, :]
    gdat.exprspec[2, :, :] = gdat.exprspec[0, :, :]

    # temp
    gdat.exprspec[where(gdat.exprspec < 0.)] = 0.

    if gdat.numbener > 1:
        gdat.exprsind = -log(gdat.exprspec[0, 1, :] / gdat.exprspec[0, 0, :]) / log(gdat.meanener[1] / gdat.meanener[0])
        gdat.exprsind[where(logical_not(isfinite(gdat.exprsind)))[0]] = 2.
    
    # temp
    gdat.exprlgal = tile(gdat.exprlgal, (3, 1)) 
    gdat.exprbgal = tile(gdat.exprbgal, (3, 1)) 
    if gdat.numbener > 1:
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

       
def retr_angldistunit(gdat, lgal, bgal, indxpixleval, retranglcosi=False):
   
    if gdat.pixltype == 'heal':
        xaxi, yaxi, zaxi = retr_unit(lgal, bgal)
        anglcosi = gdat.xaxigrid[indxpixleval] * xaxi + gdat.yaxigrid[indxpixleval] * yaxi + gdat.zaxigrid[indxpixleval] * zaxi
        
        if retranglcosi:
            return anglcosi
        else:
            angldist = arccos(anglcosi)
            return angldist
    
    else:
        angldist = sqrt((lgal - gdat.lgalgrid[indxpixleval])**2 + (bgal - gdat.bgalgrid[indxpixleval])**2)
        
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
    
    print '%5s %22s %14s %14s %14s %14s %14s %14s %10s' % ('index', 'name', 'thissamp', 'nextsamp', 'thissampvarb', 'nextsampvarb', 'diffsampvarb', 'prop', 'indxstdp')
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
        print '%5d %22s %14.4g %14.4g %14.4g %14.4g %14.4g %14s %10d' % (k, name, getattr(gdatmodi, 'thissamp')[k], getattr(gdatmodi, 'nextsamp')[k], gdatmodi.thissampvarb[k], \
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
            if listscalcomptemp[k] == 'powrslop' or listscalcomptemp[k] == 'gaus':
                if listscalcomptemp[k] == 'powrslop' or listscalcomptemp[k] == 'igam':
                    slop = gdatmodi.nextsampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[l]]
                    if listscalcomptemp[k] == 'powrslop':
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
    
    if gdat.fittnumbtrap > 0:
        gdatmodi.nextindxelemfull = deepcopy(gdatmodi.thisindxelemfull)
  
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
                gdatmodi.numbelemeval = 0
                gdatmodi.propfixp = True
            else:
                gdatmodi.numbelemeval = 2
                gdatmodi.indxtrapmodi = gdatmodi.indxsampmodi - gdat.fittindxsampcompinit
                gdatmodi.indxpoplmodi = amin(where(gdatmodi.indxtrapmodi // gdat.fittnumbtrapcumr == 0)[0])
                gdatmodi.numbparapoplinit = gdatmodi.indxtrapmodi - gdat.fittnumbtrapcuml[gdatmodi.indxpoplmodi]
                gdatmodi.indxelemmodi = [gdatmodi.numbparapoplinit // gdat.fittnumbcomp[gdatmodi.indxpoplmodi]]
                gdatmodi.indxelemfullmodi = [gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi].index(gdatmodi.indxelemmodi[0])]
                gdatmodi.indxcompmodi = gdatmodi.numbparapoplinit % gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
                gdatmodi.thisstrgcomp = gdat.fittliststrgcomp[gdatmodi.indxpoplmodi][gdatmodi.indxcompmodi]
                gdatmodi.thisliststrgcomp = [[] for l in gdat.fittindxpopl]
                gdatmodi.thisliststrgcomp[gdatmodi.indxpoplmodi].append(gdatmodi.thisstrgcomp)
                gdatmodi.thislistscalcomp = [[] for l in gdat.fittindxpopl]
                gdatmodi.thislistscalcomp[gdatmodi.indxpoplmodi].append(gdat.fittlistscalcomp[gdatmodi.indxpoplmodi][gdatmodi.indxcompmodi])
                gdatmodi.thisamplpert = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]], None]
                
            ## propose a step in a parameter
            if gdatmodi.propfixp:
                stdvpara = stdvstdp[gdat.indxstdppara[gdatmodi.indxsampmodi]]
            else:
                stdvpara = retr_propcompscal(gdat, gdatmodi, stdvstdp, gdatmodi.indxpoplmodi, gdatmodi.thisstrgcomp, indxelemfull=gdatmodi.indxelemfullmodi[0])
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, stdvpara)
            
            if gdatmodi.propfixp:
                gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_fixp(gdat, 'fitt', gdatmodi.nextsamp[gdatmodi.indxsampmodi], gdatmodi.indxsampmodi)
           
            if gdat.fittnumbtrap > 0 and gdatmodi.indxsampmodi in gdat.fittindxfixpdist:
                gdatmodi.propdist = True
                # temp
                # this should check whether rescaling is necessary
                gdatmodi.thisstrgcomp = gdat.fittnamepara[gdatmodi.indxsampmodi][-16:-12]

                ### rescale the element components due to the hyperparameter step
                rscl_elem(gdat, gdatmodi, gdatmodi.indxsampmodi)
            else:
                gdatmodi.propdist = False
            
            if not gdatmodi.propfixp:
                indxpopl = [gdatmodi.indxpoplmodi]
                
                retr_sampvarbcomp(gdat, 'fitt', gdatmodi.thisindxsampcomp, indxpopl, gdatmodi.thisliststrgcomp, gdatmodi.thislistscalcomp, gdatmodi.nextsamp, \
                                                                                                                gdatmodi.nextsampvarb, indxelemfull=gdatmodi.indxelemfullmodi[0])
        else:
            
            ## propose steps in all fixed dimensional, floating parameters
            retr_gaus(gdat, gdatmodi, gdat.indxfixpprop, stdvstdp[gdat.indxstdppara[gdat.indxfixpprop]])
            
            ### rescale the element components due to the hyperparameter steps
            if gdat.propdist:
                for k in gdat.fittindxfixpdist:
                    gdatmodi.nextsampvarb[k] = icdf_fixp(gdat, 'fitt', gdatmodi.nextsamp[k], k)
                rscl_elem(gdat, gdatmodi)

            ### element
            gdatmodi.thislfctasym = 0.
            if gdat.propcomp:
                for l in gdat.fittindxpopl:
                    for strgcomp in gdat.fittliststrgcomp[l]:
                        stdvpara = retr_propcompscal(gdat, gdatmodi, stdvstdp, l, strgcomp)
                        retr_gaus(gdat, gdatmodi, gdatmodi.thisindxsampcomp[strgcomp][l], stdvpara)

    else:
        # number of unit sample vector elements to be modified
        numbcompmodi = gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
    
    if gdatmodi.proptran:
        gdatmodi.auxipara = rand(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.indxsamptran = []
    
    if gdatmodi.propbrth or gdatmodi.propsplt:
       
        # find an empty slot in the element list
        for k in range(gdat.fittmaxmnumbpnts[gdatmodi.indxpoplmodi]):
            if not k in gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi]:
                break
       
        # sample indices to add the new element
        gdatmodi.indxsamptran.append(retr_indxsamppnts(gdat, gdatmodi.indxpoplmodi, array([k])))
        gdatmodi.indxelemmodi = [k]
        gdatmodi.nextindxelemfull[gdatmodi.indxpoplmodi].append(k)
        gdatmodi.indxelemfullmodi = [gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]].astype(int)]
    if gdatmodi.propbrth:
        
        gdatmodi.numbelemeval = 1
        
        # sample auxiliary variables
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.auxipara
    
    # death
    if gdatmodi.propdeth:
        
        gdatmodi.numbelemeval = 1
        
        # occupied element index to be killed
        dethindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # element index to be killed
        gdatmodi.indxelemmodi = []
        gdatmodi.indxelemfullmodi = []
        gdatmodi.indxelemmodi.append(gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][dethindxindxpnts])
        gdatmodi.indxelemfullmodi.append(dethindxindxpnts)
        # sample indices to add the new element
        gdatmodi.indxsamptran.append(gdat.fittindxsampcompinit + gdat.fittnumbtrapcuml[gdatmodi.indxpoplmodi] + \
                                                    gdatmodi.indxelemmodi[0] * gdat.fittnumbcomp[gdatmodi.indxpoplmodi] + gdat.fittindxcomp[gdatmodi.indxpoplmodi])
    
    if gdatmodi.propsplt or gdatmodi.propmerg:
        gdatmodi.comppare = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.compfrst = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.compseco = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
    
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbelemeval = 3
        
        gdatmodi.indxelemfullsplt = choice(arange(gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        gdatmodi.indxelemsplt = gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][gdatmodi.indxelemfullsplt]
        gdatmodi.indxelemfullmodi.insert(0, gdatmodi.indxelemfullsplt)
        gdatmodi.indxelemmodi.insert(0, gdatmodi.indxelemsplt)

        # boundaries of the sample vector indices to be modified
        ## first
        gdatmodi.indxsampfrst = gdat.fittindxsampcompinit + gdat.fittnumbtrap * gdatmodi.indxpoplmodi + gdatmodi.indxelemsplt * gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        gdatmodi.indxsampfrstfinl = gdatmodi.indxsampfrst + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        gdatmodi.indxsamptran.insert(0, arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl))
        
        ## second
        gdatmodi.indxsampseco = gdat.fittindxsampcompinit + gdat.fittnumbtrap * gdatmodi.indxpoplmodi + gdatmodi.indxelemmodi[1] * gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        gdatmodi.indxsampsecofinl = gdatmodi.indxsampseco + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
        
        # concatenated sample vector indices to be modified
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl, dtype=int), \
                                                                                    arange(gdatmodi.indxsampseco, gdatmodi.indxsampsecofinl, dtype=int)))
        
        for strgcomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
            comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
            setattr(gdatmodi, strgcomp + 'pare', comp)
        gdatmodi.comppare[2] = getattr(gdatmodi, gdat.namecompampl + 'pare')

        # determine the new element parameters
        # temp -- only valid for power-law energy spectrum
        for k, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
            if k == 0 or k == 1:
                gdatmodi.auxipara[k] = randn() * gdat.radispmr
            elif k == 2:
                gdatmodi.auxipara[k] = rand()
            else:
                if gdat.fittlistscalcomp[gdatmodi.indxpoplmodi] == 'self':
                    minm = getattr(gdat, 'fittminm' + strgcomp)
                    fact = getattr(gdat, 'fittfact' + strgcomp)
                    gdatmodi.auxipara[k] = icdf_self(rand(), minm, fact)
                elif gdat.fittlistscalcomp[gdatmodi.indxpoplmodi] == 'gaus':
                    distmean = sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[gdatmodi.indxpoplmodi]]
                    diststdv = sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[gdatmodi.indxpoplmodi]]
                    gdatmodi.auxipara[k] = icdf_gaus(rand(), distmean, diststdv)
        
        gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.auxipara[2]) * gdatmodi.auxipara[0]
        gdatmodi.compfrst[1] = gdatmodi.comppare[1] + (1. - gdatmodi.auxipara[2]) * gdatmodi.auxipara[1]
        gdatmodi.compfrst[2] = gdatmodi.auxipara[2] * gdatmodi.comppare[2]
        gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.auxipara[2] * gdatmodi.auxipara[0]
        gdatmodi.compseco[1] = gdatmodi.comppare[1] - gdatmodi.auxipara[2] * gdatmodi.auxipara[1]
        gdatmodi.compseco[2] = (1. - gdatmodi.auxipara[2]) * gdatmodi.comppare[2]
        #for k in range(gdat.fittnumbcomp[gdatmodi.indxpoplmodi]):
        #    if k > 2:
        #        gdatmodi.compfrst[k] = gdatmodi.comppare[k]
        #        gdatmodi.compseco[k] = gdatmodi.comppare[k]
        
        
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0][:3]] = cdfn_comp(gdat, gdatmodi, gdatmodi.compfrst)[:3]
        gdatmodi.nextsamp[gdatmodi.indxsamptran[1][:3]] = cdfn_comp(gdat, gdatmodi, gdatmodi.compseco)[:3]

        if fabs(gdatmodi.compfrst[0]) > gdat.maxmgang or fabs(gdatmodi.compseco[0]) > gdat.maxmgang or \
                                      fabs(gdatmodi.compfrst[1]) > gdat.maxmgang or fabs(gdatmodi.compseco[1]) > gdat.maxmgang or \
                                      gdatmodi.compfrst[2] < getattr(gdat, 'fittminm' + gdat.namecompampl) or gdatmodi.compseco[2] < getattr(gdat, 'fittminm' + gdat.namecompampl):
            if gdat.verbtype > 1:
                print 'Proposal rejected due to component falling outside the prior.'
                print 'gdat.maxmgang'
                print gdat.maxmgang * gdat.anglfact
                print 'getattr(gdat, fittminm + gdat.namecompampl)'
                print getattr(gdat, 'fittminm' + gdat.namecompampl)
                print 'gdatmodi.compfrst[2]'
                print gdatmodi.compfrst[2]
                print 'gdatmodi.compseco[2]'
                print gdatmodi.compseco[2]
                print

            gdatmodi.thisaccpprio = False

        # calculate the list of pairs
        if gdatmodi.thisaccpprio:

            # calculate the list of pairs
            ## proposed
            lgal = concatenate((array([gdatmodi.compfrst[0], gdatmodi.compseco[0]]), \
                                                           setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], gdatmodi.comppare[0])))
            bgal = concatenate((array([gdatmodi.compfrst[1], gdatmodi.compseco[1]]), \
                                                           setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]], gdatmodi.comppare[1])))
            
            if gdat.verbtype > 1:
                print 'Element positions prior to calculating the pair list for the reverse move of a split'
                print 'lgal'
                print lgal * gdat.anglfact
                print 'bgal'
                print bgal * gdat.anglfact

            #gdatmodi.thislistpair = retr_listpair(gdat, lgal, bgal)
            gdatmodi.thislistpair = [0]
            gdatmodi.thisnumbpair = len(gdatmodi.thislistpair)

            if gdatmodi.thisnumbpair == 0:
                raise Exception('Number of pairs should not be zero in the reverse proposal of a split')

    if gdatmodi.propmerg:
        
        # number of point sources to be modified
        gdatmodi.numbelemeval = 3
        
        # calculate the current list of pairs
        # temp
        #gdatmodi.thislistpair = retr_listpair(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], \
        #                                                                                gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]])
        gdatmodi.thislistpair = [0]
        gdatmodi.thisnumbpair = len(gdatmodi.thislistpair)
         
        if gdat.verbtype > 1:
            print 'thislistpair'
            print gdatmodi.thislistpair
           
        # check if merge will be proposed
        if gdatmodi.thisnumbpair == 0:
            gdatmodi.thisaccpprio = False
            if gdat.verbtype > 1:
                print 'Proposal rejected due to not finding a pair to merge.'
        else:

            # sample a pair
            #indxpairtemp = choice(arange(gdatmodi.thisnumbpair))
    
            # determine element indices to be merged
            #gdatmodi.indxelemfullmergfrst = gdatmodi.thislistpair[indxpairtemp][0]
            #gdatmodi.indxelemfullmergseco = gdatmodi.thislistpair[indxpairtemp][1]
            #gdatmodi.indxelemfullmodi = [gdatmodi.indxelemfullmergfrst, gdatmodi.indxelemfullmergseco]
            
            gdatmodi.indxelemfullmodi = sort(choice(arange(len(gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi])), size=2, replace=False))
            gdatmodi.indxelemfullmergfrst = gdatmodi.indxelemfullmodi[0]
            gdatmodi.indxelemfullmergseco = gdatmodi.indxelemfullmodi[1]

            ## first element index to be merged
            gdatmodi.mergindxpntsfrst = gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmergfrst]

            ## second element index to be merged
            gdatmodi.mergindxpntsseco = gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmergseco]
            
            # determine indices of the modified elements in the sample vector
            ## first element
            # temp -- this would not work for multiple populations !
            gdatmodi.indxsampfrst = gdat.fittindxsampcompinit + gdat.fittnumbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxpntsfrst
            gdatmodi.indxsampfrstfinl = gdatmodi.indxsampfrst + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
            gdatmodi.indxsamptran.append(arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl))

            ## second element
            gdatmodi.indxsampseco = gdat.fittindxsampcompinit + gdat.fittnumbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxpntsseco
            gdatmodi.indxsampsecofinl = gdatmodi.indxsampseco + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
            gdatmodi.indxsamptran.append(arange(gdatmodi.indxsampseco, gdatmodi.indxsampsecofinl))
            
            # indices of the sample vector elements to be modified
            gdatmodi.indxsampmodi = arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl)

            # indices of the element to be merged
            gdatmodi.indxelemmodi = [gdatmodi.mergindxpntsfrst, gdatmodi.mergindxpntsseco]

            # parameters of the elements to be merged
            for k, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
                ## first
                gdatmodi.compfrst[k] = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmergfrst]]
                ## second
                gdatmodi.compseco[k] = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmergseco]]

            # auxiliary parameters
            gdatmodi.auxipara[0] = gdatmodi.compseco[0] - gdatmodi.compfrst[0]
            gdatmodi.auxipara[1] = gdatmodi.compseco[1] - gdatmodi.compfrst[1]
            gdatmodi.auxipara[2] = gdatmodi.compfrst[2] / (gdatmodi.compfrst[2] + gdatmodi.compseco[2]) 
            for k, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
                if k > 2:
                    gdatmodi.auxipara[k] = gdatmodi.compseco[k]

            # merged element
            gdatmodi.comppare[2] = gdatmodi.compfrst[2] + gdatmodi.compseco[2]
            if gdatmodi.comppare[2] > getattr(gdat, 'maxm' + gdat.namecompampl):
                gdatmodi.thisaccpprio = False
                if gdat.verbtype > 1:
                    print 'Proposal rejected due to falling outside the prior.'
            gdatmodi.comppare[0] = gdatmodi.compfrst[0] + (1. - gdatmodi.auxipara[2]) * (gdatmodi.compseco[0] - gdatmodi.compfrst[0])
            gdatmodi.comppare[1] = gdatmodi.compfrst[1] + (1. - gdatmodi.auxipara[2]) * (gdatmodi.compseco[1] - gdatmodi.compfrst[1])
            for k, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]):
                if k < 2:
                    gdatmodi.comppare[k] = gdatmodi.compfrst[k] + (1. - gdatmodi.auxipara[2]) * (gdatmodi.compseco[k] - gdatmodi.compfrst[k])
                elif k == 2:
                    gdatmodi.comppare[k] = gdatmodi.compfrst[k] + gdatmodi.compseco[k]
                else:
                    setattr(gdatmodi, strgcomp + 'pare', strgcomp + 'frst')

            gdatmodi.nextsamp[gdatmodi.indxsamptran[0][:3]] = cdfn_comp(gdat, gdatmodi, gdatmodi.comppare)[:3]
            gdatmodi.nextsamp[gdatmodi.indxsamptran[0][:3]] = gdatmodi.thissamp[gdatmodi.indxsamptran[0][3:]]
            gdatmodi.nextsampvarb[gdatmodi.indxsamptran[0]] = gdatmodi.comppare

            # calculate the proposed list of pairs
            if gdat.verbtype > 1:
                print 'mergindxfrst: ', gdatmodi.mergindxpntsfrst
                print 'gdatmodi.indxelemfullmergfrst: ', gdatmodi.indxelemfullmergfrst
                print 'mergindxseco: ', gdatmodi.mergindxpntsseco
                print 'gdatmodi.indxelemfullmergseco: ', gdatmodi.indxelemfullmergseco
                print 'indxsampfrst: ', gdatmodi.indxsampfrst
                print 'gdatmodi.indxsampfrstfinl: ', gdatmodi.indxsampfrstfinl
                print 'indxsampseco: ', gdatmodi.indxsampseco
                print 'gdatmodi.indxsampsecofinl: ', gdatmodi.indxsampsecofinl
            
            #lgal = concatenate((array([gdatmodi.comppare[0]]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi]], \
            #                                                                                            array([gdatmodi.compfrst[0], gdatmodi.compseco[0]]))))
            #bgal = concatenate((array([gdatmodi.comppare[1]]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi]], \
            #                                                                                            array([gdatmodi.compfrst[1], gdatmodi.compseco[1]]))))
        
    if gdat.verbtype > 1 and (gdatmodi.propsplt or gdatmodi.thisaccpprio and gdatmodi.propmerg):
        print 'lgalfrst: ', gdat.anglfact * gdatmodi.compfrst[0]
        print 'bgalfrst: ', gdat.anglfact * gdatmodi.compfrst[1]
        print 'amplfrst: ', gdatmodi.compfrst[2]
        print 'lgalseco: ', gdat.anglfact * gdatmodi.compseco[0]
        print 'bgalseco: ', gdat.anglfact * gdatmodi.compseco[1]
        print 'amplseco: ', gdatmodi.compseco[2]
        print 'lgalpare: ', gdat.anglfact * gdatmodi.comppare[0]
        print 'bgalpare: ', gdat.anglfact * gdatmodi.comppare[1]
        print 'fluxpare: ', gdatmodi.comppare[2]
        print 'auxipara[0]: ', gdat.anglfact * gdatmodi.auxipara[0]
        print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
        print 'auxipara[2]: ', gdatmodi.auxipara[2]
                
    # change the number of elements
    if gdatmodi.propbrth or gdatmodi.propsplt:
        gdatmodi.nextsamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] + 1
    if gdatmodi.propdeth or gdatmodi.propmerg:
        gdatmodi.nextsamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]] - 1
    
    if (gdatmodi.propdeth or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        # remove the element from the occupied element list
        for a, indxpnts in enumerate(gdatmodi.indxelemmodi):
            if a == 0 and gdatmodi.propdeth or a == 1 and gdatmodi.propmerg:
                gdatmodi.nextindxelemfull[gdatmodi.indxpoplmodi].remove(indxpnts)
    
    if gdat.fittnumbtrap > 0:
        if gdatmodi.propwith:
            if gdat.propwithsing:
                if gdatmodi.propdist:
                    gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampmodi]), gdatmodi.thisindxsampcomp[gdatmodi.thisstrgcomp][gdatmodi.indxpoplmodi]))
            else:
                gdatmodi.indxsampmodi = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp['comp'])))
            gdatmodi.indxsampchec = gdatmodi.indxsampmodi
        else:
            gdatmodi.indxsampchec = []
            if (gdatmodi.propbrth or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[0]))
            if gdatmodi.propdeth:
                gdatmodi.indxsampmodi = gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None]
            if gdatmodi.propsplt:
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[0], gdatmodi.indxsamptran[1]))
    else:
        gdatmodi.indxsampchec = gdat.indxfixpprop
        gdatmodi.indxsampmodi = gdat.indxfixpprop
    
    # reject the sample if proposal is outside the prior
    indxchecfail = where((gdatmodi.nextsamp[gdatmodi.indxsampchec] < 0.) | (gdatmodi.nextsamp[gdatmodi.indxsampchec] > 1.))[0]
    if indxchecfail.size > 0:
        if gdat.verbtype > 1:
            for k in range(20):
                print 'Proposal rejected due to proposal outside the prior during the common check'
            print 'indxchecfail'
            print indxchecfail
            print
        gdatmodi.thisaccpprio = False
    
    if gdat.verbtype > 1:
        print 'gdatmodi.thisindxelemfull'
        print gdatmodi.thisindxelemfull
        print 'gdatmodi.nextindxelemfull'
        print gdatmodi.nextindxelemfull
        print 'gdatmodi.thisaccpprio'
        print gdatmodi.thisaccpprio
        print 'gdatmodi.indxsampchec'
        print gdatmodi.indxsampchec
        if not gdatmodi.propfixp and not (gdatmodi.propmerg and not gdatmodi.thisaccpprio):
            print 'gdatmodi.indxelemmodi'
            print gdatmodi.indxelemmodi
            print 'gdatmodi.indxelemfullmodi'
            print gdatmodi.indxelemfullmodi
        if gdatmodi.proptran:
            print 'gdatmodi.indxsamptran'
            for indxsamptran in gdatmodi.indxsamptran:
                print indxsamptran

    if gdatmodi.propwith:
        gdatmodi.nextindxsampcomp = gdatmodi.thisindxsampcomp
    else:
        gdatmodi.nextindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.nextindxelemfull, 'fitt')
    
    # temp -- this is inefficient for propwithsing proposals
    if gdatmodi.thisaccpprio and (gdatmodi.propbrth or gdatmodi.propsplt or gdatmodi.propwith and not gdat.propwithsing):
        gdatmodi.nextsampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.nextsamp, gdatmodi.nextindxsampcomp)
   
    if gdat.verbtype > 1:
        show_samp(gdat, gdatmodi)
    
    if gdat.diagmode:
        for strgcomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
            if gdatmodi.propbrth and (gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi]] != \
                                            gdatmodi.nextsampvarb[gdatmodi.thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi]]).any():
                raise Exception('')

    if gdat.propwithsing and not gdatmodi.propfixp:
        gdatmodi.dicttemp = dict()
        for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
            if gdatmodi.propwith:
                thiscomp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                nextcomp = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                if namecomp == gdat.namefeateval:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([-thiscomp, nextcomp])
                else:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([thiscomp, nextcomp])
            elif gdatmodi.propbrth:
                comp = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                if namecomp == gdat.namefeateval:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([comp])
                else:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([comp])
            elif gdatmodi.propdeth:
                comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                if namecomp == gdat.namefeateval:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([-comp])
                else:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([comp])
            elif gdatmodi.propsplt:
                comppare = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                compfrst = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                compseco = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[1]]]
                if namecomp == gdat.namefeateval:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([-comppare, compfrst, compseco])
                else:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([comppare, compfrst, compseco])
            elif gdatmodi.propmerg and gdatmodi.thisaccpprio:
                compfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                compseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[1]]]
                comppare = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmodi[0]]]
                if namecomp == gdat.namefeateval:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([-compfrst, -compseco, comppare])
                else:
                    gdatmodi.dicttemp[namecomp + 'eval'] = array([compfrst, compseco, comppare])
        
        if gdat.verbtype > 1:
            print 'gdatmodi.dicttemp'
            print gdatmodi.dicttemp
            print

    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        
        ## Jacobian
        if gdatmodi.propsplt:
            gdatmodi.thisljcbfact = log(gdatmodi.comppare[2])
        else:
            gdatmodi.thisljcbfact = -log(gdatmodi.comppare[2])
         
        ## combinatorial factor
        # temp
        #thisnumbpnts = gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]]
        #if gdatmodi.propsplt:
        #    gdatmodi.thislcomfact = log(gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2 / gdatmodi.thisnumbpair)
        #else:
        #    gdatmodi.thislcomfact = log(gdatmodi.thisnumbpair / gdatmodi.thissamp[gdat.fittindxfixpnumbpnts[gdatmodi.indxpoplmodi]]**2)
        gdatmodi.thislcomfact = 0.


def retr_propcompscal(gdat, gdatmodi, stdvstdp, l, strgcomp, indxelemfull=slice(None)):
    
    thiscompampl = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
    nextcompampl = gdatmodi.nextsampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
    minmcompampl = getattr(gdat, 'minm' + gdat.namecompampl)
    stdvstdpcomp = stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)]
    thiscompunit = gdatmodi.thissamp[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
    nextcompunit = gdatmodi.nextsamp[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
    if strgcomp == gdat.namecompampl:
        # temp -- this only works if compampl is powr distributed
        gdatmodi.thisstdv = stdvstdpcomp / (thiscompampl / minmcompampl)**2.
        gdatmodi.nextstdv = stdvstdpcomp / (nextcompampl / minmcompampl)**2.
        gdatmodi.thislfctasym += sum(0.5 * (nextcompunit - thiscompunit)**2 * (1. / gdatmodi.thisstdv**2 - 1. / gdatmodi.nextstdv**2))
    else:
        gdatmodi.thisstdv = stdvstdpcomp / (minimum(thiscompampl, nextcompampl) / minmcompampl)**0.5

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
        print 'gdat.listindxelemfull'
        print gdat.listindxelemfull

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
            for k in range(len(gdat.listindxelemfull[n][l])):
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
        indxsampcomp = retr_indxsampcomp(gdat, gdat.listindxelemfull[n], 'fitt') 
        for l in gdat.fittindxpopl:
            for k in arange(len(gdat.listindxelemfull[n][l])):
                for m, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                    arrystks[indxstks[n][l][k], m] = gdat.listsamp[n, indxsampcomp[strgcomp][l][k]]

    if gdat.verbtype > 0:
        print 'Constructing the distance matrix for %d stacked samples...' % arrystks.shape[0]
        timeinit = gdat.functime()
    
    gdat.distthrs = empty(gdat.fittnumbcomptotl)
    for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
        gdat.distthrs[k] = gdat.stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)]
    
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
            dist = fabs(arrystks[n, k] - arrystks[:, k])
            indxstks = where(dist < gdat.distthrs[k])[0]
            if indxstks.size > 0:
                for j in indxstks:
                    cntr += 1
                    listdisttemp[k].append(dist[j])
                    indxstksrows[k].append(n)
                    indxstkscols[k].append(j)
            
            nextperc = floor(100. * float(k * numbstks + n) / numbstks / gdat.fittnumbcomptotl)
            if nextperc > thisperc:
                thisperc = nextperc
                print '%3d%% completed.' % thisperc
            if cntr > 1e6:
                break
        
        listdisttemp[k] = array(listdisttemp[k])
        indxstksrows[k] = array(indxstksrows[k])
        indxstkscols[k] = array(indxstkscols[k])

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
        numbdist = zeros(numbstks, dtype=int) - 1
        for p in range(len(indxstksleft)):
            indxindx = where((listdist[0][indxstksleft[p], :].toarray().flatten() * 2. * gdat.maxmlgal < gdat.anglassc) & \
                             (listdist[1][indxstksleft[p], :].toarray().flatten() * 2. * gdat.maxmbgal < gdat.anglassc))[0]
            numbdist[indxstksleft[p]] = indxindx.size
            
        prvlmaxmesti = amax(numbdist) / float(gdat.numbsamptotl)
        
        if False and gdat.verbtype > 0 and cntr % 10 == 0:
            print 'cntr'
            print cntr
            print 'numbdist'
            print numbdist
            print 'gdat.numbsamptotl'
            print gdat.numbsamptotl

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
                #    gdatmodi.thisindxelemfull = deepcopy(gdat.listindxelemfull[n])
                #    for r in range(len(indxstksassc)): 
                #        calc_poststkscond(gdat, indxstksassc)
                #    gdatmodi.thisindxelemfull = [[] for l in gdat.fittindxpopl]
                #    for indxstkstemp in indxstksleft:
                #        indxsamptotlcntr = indxtupl[indxstkstemp][0]
                #        indxpoplcntr = indxtupl[indxstkstemp][1]
                #        indxpntscntr = indxtupl[indxstkstemp][2]
                #        gdatmodi.thissampvarb = gdat.listsampvarb[indxsamptotlcntr, :]
                #        gdatmodi.thisindxelemfull[].append()

                #    plot_genemaps(gdat, gdatmodi, 'this', 'datacnts', thisindxener=0, thisindxevtt=0, cond=True)
                
                if gdat.verbtype > 1:
                    print 

            cntr += 1
        
        if gdat.verbtype > 1:
            print 
            print 
            print 
        
    if gdat.verbtype > 1:
        print 'gdat.listindxelemfull'
        print gdat.listindxelemfull
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
    gdat.minmprvl = 0.
    gdat.maxmprvl = 1.
    retr_axis(gdat, 'prvl', gdat.minmprvl, gdat.maxmprvl, 10)
    gdat.histprvl = histogram(gdat.prvl, bins=gdat.binsprvl)[0]
    if gdat.makeplot:
        pathcond = getattr(gdat, 'path' + gdat.namesampdist + 'finlcond')
        for k, strgcomp in enumerate(gdat.fittliststrgcomptotl):
            path = pathcond + 'histdist' + strgcomp 
            listtemp = copy(listdist[k].toarray()).flatten()
            listtemp = listtemp[where(listtemp != 1e20)[0]]
            tdpy.mcmc.plot_hist(path, listtemp, r'$\Delta \tilde{' + getattr(gdat, 'labl' + strgcomp) + '}$')
            path = pathcond + 'histprvl'
            tdpy.mcmc.plot_hist(path, gdat.prvl, r'$p$')
    gdat.prvlthrs = 0.1 
    gdat.indxprvlhigh = where(gdat.prvl > gdat.prvlthrs)[0]
    gdat.numbprvlhigh = gdat.indxprvlhigh.size


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
    
    factmcutfromdefs = pi * adishost**2 * critmden * asca * retr_mcutfrommscl(fracacutasca)

    return factmcutfromdefs


def retr_mcut(gdat, defs, asca, acut):
    
    mscl = defs * pi * gdat.adishost**2 * gdat.critmden * asca
    fracacutasca = acut / asca
    mcut = mscl * retr_mcutfrommscl(fracacutasca)
    
    return mcut


def retr_mcutfrommscl(fracacutasca):
    
    mcut = fracacutasca**2 / (fracacutasca**2 + 1.)**2 * ((fracacutasca**2 - 1.) * log(fracacutasca) + fracacutasca * pi - (fracacutasca**2 + 1.))

    return mcut


def retr_negalogt(varb):
    
    negalogt = sign(varb) * log10(fabs(varb))
    
    return negalogt


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
    
    gdat.liststrgmodl = ['fitt']
    if gdat.datatype == 'mock':
        gdat.liststrgmodl += ['true']
    
    gdat.listnamefeateval = ['lgal', 'bgal', 'spec']
        
    if gdat.elemtype == 'lght':
        gdat.liststrgfeatplot = []
    if gdat.elemtype == 'lens':
        #gdat.liststrgfeatplot = ['mcutcorr']
        gdat.liststrgfeatplot = []

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
            
        for nameseco in ['histodim', 'histtdim', 'assc', 'scattdim', 'cmpl', 'fdis']:
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
        gdat.ascaglob = 0.05 / gdat.anglfact
        gdat.acutglob = 1. / gdat.anglfact
    
    gdat.cutfdefs = 3e-3 / gdat.anglfact

    if gdat.elemtype == 'lght':
        gdat.namefeateval = 'flux'
    if gdat.elemtype == 'lens':
        gdat.namefeateval = 'defs'
    
    # plotting
    gdat.legdsampdist = 'Posterior'
    gdat.legdsamp = 'Sample'
    gdat.legdmlik = 'Maximum likelihood'
    gdat.legdmedi = 'Median'
    gdat.legdstdv = 'Std. dev.'
    
    # optimization period
    gdat.numbswepoptiprop = 10 * gdat.fittnumbpara

    # p value contours 
    gdat.pvalcont = [0.317, 0.0455, 2.7e-3, 6e-5, 1.3e-6]

    ## number of bins in histogram plots
    gdat.numbbinsplot = 20
    gdat.indxbinsplot = arange(gdat.numbbinsplot)
    
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

        gdat.lablfracsubhdeltbein = '$f_{E,sub}$'
        gdat.scalfracsubhdeltbein = 'logt'
        
        gdat.lablfracsubhintgbein = '$f_{<E,sub}$'
        gdat.scalfracsubhintgbein = 'logt'
        
        gdat.lablmasssubhdeltbein = '$M_{E,sub}$'
        gdat.scalmasssubhdeltbein = 'logt'
        
        gdat.lablmasssubhintgbein = '$M_{<E,sub}$'
        gdat.scalmasssubhintgbein = 'logt'
        
        gdat.lablmasshostdeltbein = '$M_{E,hst}$'
        gdat.scalmasshostdeltbein = 'logt'
    
        gdat.lablmasshostintgbein = '$M_{<E,hst}$'
        gdat.scalmasshostintgbein = 'logt'
    
    gdat.scalmaxmnumbpnts = 'logt'
    gdat.scalmedilliktotl = 'logt'

    gdat.lablener = 'E'
    #gdat.lablenertotl = '$%s$ [%s]' % (gdat.lablener, gdat.strgenerunit)
    
    gdat.lablgang = r'\theta'
    gdat.lablaang = r'\phi'
    gdat.labllgalunit = gdat.lablgangunit
    gdat.lablbgalunit = gdat.lablgangunit
   
    gdat.lablanglfromhost = r'\theta_{0,hst}'
    gdat.lablanglfromhostunit = gdat.lablgangunit
    gdat.lablmassintg = r'<M>_{<r}'
    gdat.lablmassintgunit = r'$M_{\odot}$'
    gdat.lablmassdelt = r'<M>_r'
    gdat.lablmassdeltunit = r'$M_{\odot}$'
    gdat.lablfracsubhintg = r'<f>_{<r,sub}'
    gdat.lablfracsubhdelt = r'<f>_{r,sub}'

    gdat.labldefs = r'\alpha_s'
    gdat.lablflux = 'f'
    gdat.lablnobj = 'p'
    
    gdat.labldeflprof = r'\alpha_a'
    gdat.labldeflprofunit = u'$^{\prime\prime}$'
    
    gdat.labldefsunit = u'$^{\prime\prime}$'
    if gdat.exprtype == 'ferm':
        gdat.lablfluxunit = 'cm$^{-2}$ s$^{-1}$ GeV$^{-1}$'
        gdat.lablfluxsoldunit = 'cm$^{-2}$ s$^{-1}$ GeV$^{-1} sr$^{-1}$'
    if gdat.exprtype == 'chan':
        gdat.lablfluxunit = 'cm$^{-2}$ s$^{-1}$ KeV$^{-1}$'
        gdat.lablfluxsoldunit = 'cm$^{-2}$ s$^{-1}$ KeV$^{-1} sr$^{-1}$'
    if gdat.exprtype == 'hubb':
        gdat.lablfluxunit = r'erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$'
        gdat.lablfluxsoldunit = r'erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$ sr$^{-1}$'
    
    gdat.lablfluxsold = 'f'
    
    for l in gdat.trueindxpopl:
        setattr(gdat, 'lablcmplpop%d' % l, '$c_{%d}$' % l)
        setattr(gdat, 'lablfdispop%d' % l, '$f_{%d}$' % l)

    gdat.lablprvl = '$p$'
    
    gdat.lablsind = 's'
    gdat.lablcurv = r'\kappa'
    gdat.lablexpo = r'\epsilon'
    gdat.lablexpounit = gdat.strgenerunit
    
    gdat.labllliktotl = r'\mathcal{L}'
    gdat.lablmedilliktotl = r'\mathcal{L}^{med}'
    
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
    gdat.labldeltllik = r'\Delta_n \ln \mathcal{L}'
    gdat.labldiss = r'\theta_{sa}'
    gdat.labldissunit = gdat.lablgangunit
    
    gdat.lablrele = r'\langle|\vec{\alpha}_n \cdot \vec{\nabla} k_m| \rangle'
    
    gdat.lablrelc = r'\langle\vec{\alpha}_n \cdot \vec{\nabla} k_m \rangle'
    
    gdat.lablreld = r'\langle|\vec{\alpha}_n \cdot \vec{\nabla} k_d| \rangle'
    
    gdat.lablreln = r'\langle \Delta \theta_{pix} |\hat{\alpha}_n \cdot \vec{\nabla} k_m| / \alpha_{s,n} \rangle'
    
    # define the labels for the selection features
    for namefeat in gdat.listnamefeatsele:
        for namesele in gdat.listnamesele:
            setattr(gdat, 'labl' + namefeat + namesele, getattr(gdat, 'labl' + namefeat))
            # temp
            #setattr(gdat, 'labl' + namefeat + namesele + 'unit', getattr(gdat, 'labl' + namefeat + 'unit'))
    
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
                setattr(gdat, 'labl' + name + 'unit', '')
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
    
    gdat.minmexpo = 1.
    gdat.maxmexpo = 10.
    
    for l in gdat.trueindxpopl:
        setattr(gdat, 'minmcmplpop%d' % l, 0.)
        setattr(gdat, 'maxmcmplpop%d' % l, 1.)
        setattr(gdat, 'scalcmplpop%d' % l, 'self')

    gdat.minmdeltllik = 1.
    gdat.maxmdeltllik = 1e3
    gdat.minmdiss = 0.
    gdat.maxmdiss = gdat.maxmgang * sqrt(2.)
    
    gdat.minmrele = 0.001
    gdat.maxmrele = 1.

    gdat.minmreln = 0.3
    gdat.maxmreln = 10.

    gdat.minmreld = 0.1
    gdat.maxmreld = 10.

    gdat.minmrelc = 0.01
    gdat.maxmrelc = 1.

    gdat.minmmcut = 1e7
    gdat.maxmmcut = 5e9
    gdat.minmmcutcorr = gdat.minmmcut
    gdat.maxmmcutcorr = gdat.maxmmcut

    gdat.minmbein = 0.
    gdat.maxmbein = 1. / gdat.anglfact
    
    # scalar variables
    gdat.minmdeflprof = 1e-3 / gdat.anglfact
    gdat.maxmdeflprof = 0.1 / gdat.anglfact

    gdat.minmfracsubh = 0.
    gdat.maxmfracsubh = 0.3
    gdat.scalfracsubh = 'self'
    # temp -- automize this
    gdat.minmfracsubhcorr = gdat.minmfracsubh
    gdat.maxmfracsubhcorr = gdat.maxmfracsubh
    gdat.scalfracsubhcorr = gdat.scalfracsubh
    
    gdat.minmmasshost = 1e10
    gdat.maxmmasshost = 1e13
    gdat.scalmasshost = 'self'
    
    gdat.minmmasssubh = 1e8
    gdat.maxmmasssubh = 1e10
    gdat.scalmasssubh = 'self'

    # set up the indices of the fitting model
    retr_indxsamp(gdat)
    
    if gdat.datatype == 'inpt':
        for l in gdat.fittindxpopl:
            for strgpdfn in gdat.fittliststrgpdfnprio[l]:
                if strgpdfn.startswith('gaum') and gdat.fittlgalprio == None and gdat.fittbgalprio == None:
                    raise Exception('If spatdisttype is "gaus", spatial coordinates of the prior catalog should be provided via lgalprio and bgalprio.')
    # temp
    for strgmodl in gdat.liststrgmodl:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
        for l in indxpopl:
            for strgfeat, strgpdfn in zip(liststrgfeatprio[l], liststrgpdfnprio[l]):
                if strgpdfn == 'self':
                    minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                    maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                    fact = maxm - minm
                    setattr(gdat, strgmodl + 'fact' + strgfeat, fact)
    
    if gdat.propcomp == None:
        gdat.propcomp = gdat.fittnumbtrap > 0

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
    if gdat.elemtype == 'lens':
        for strgtemp in ['delt', 'intg']:
            gdat.listnamevarbscal += ['masshost' + strgtemp + 'bein']
            setattr(gdat, 'minmmasshost' + strgtemp + 'bein', gdat.minmmasshost)
            setattr(gdat, 'maxmmasshost' + strgtemp + 'bein', gdat.maxmmasshost)
    if gdat.fittnumbtrap > 0:
        for l in gdat.trueindxpopl:
            gdat.listnamevarbscal += ['cmplpop%d' % l]
        for l in gdat.fittindxpopl:
            gdat.listnamevarbscal += ['fdispop%d' % l]
        if gdat.elemtype == 'lens':
            for strgtemp in ['delt', 'intg']:
                gdat.listnamevarbscal += ['masssubh' + strgtemp + 'bein', 'fracsubh' + strgtemp + 'bein'] 
                setattr(gdat, 'minmmasssubh' + strgtemp + 'bein', gdat.minmmasssubh)
                setattr(gdat, 'maxmmasssubh' + strgtemp + 'bein', gdat.maxmmasssubh)
                setattr(gdat, 'minmfracsubh' + strgtemp + 'bein', gdat.minmfracsubh)
                setattr(gdat, 'maxmfracsubh' + strgtemp + 'bein', gdat.maxmfracsubh)
            if gdat.priofactdoff != 0. or (gdat.fittnamefixp[gdat.fittindxfixpmeanpnts] == 'logt').any():
                gdat.listnamevarbscal += ['fracsubhcorr']
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
        for strgfeat in liststrgfeattotl + gdat.liststrgfeatplot:
            labl = getattr(gdat, 'labl' + strgfeat)
            lablunit = getattr(gdat, 'labl' + strgfeat + 'unit')
            if lablunit != '':
                setattr(gdat, 'labl' + strgfeat + 'unit', ' [%s]' % lablunit)
            
            labltotl = '$%s$%s' % (labl, lablunit)
            setattr(gdat, 'labl' + strgfeat + 'totl', labltotl)
            if strgfeat.startswith('defs') or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal' or \
                                                                                                        strgfeat == 'diss' or strgfeat == 'asca' or strgfeat == 'acut':
                setattr(gdat, 'fact' + strgfeat + 'plot', gdat.anglfact)
            else:
                setattr(gdat, 'fact' + strgfeat + 'plot', 1.)
            
            if strgfeat.startswith(gdat.namefeatsign) or strgfeat.startswith(gdat.namecompampl) or strgfeat == 'expo' or strgfeat == 'cnts' or \
                                    strgfeat == 'relc' or strgfeat == 'rele' or strgfeat == 'reln' or strgfeat == 'reld' or strgfeat.startswith('mcut') or strgfeat == 'deltllik':
                setattr(gdat, 'scal' + strgfeat + 'plot', 'logt')
            else:
                setattr(gdat, 'scal' + strgfeat + 'plot', 'self')
    
    # log-prior register
    ## indices of penalization term
    indxlpripena = 0
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    numb = 0
    for l in gdat.fittindxpopl:
        numb += len(gdat.fittliststrgfeat[l])
    if gdat.fittmaxmnumbpntstotl > 0:
        gdat.numblpri = 1 + gdat.fittnumbpopl + numb
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
    gdat.numbangl = 1000
    gdat.numbanglhalf = 10
    gdat.indxanglhalf = arange(gdat.numbanglhalf)
    gdat.numbanglfull = 1000
    gdat.indxanglfull = arange(gdat.numbanglfull)
    if gdat.evalcirc == 'full':
        gdat.maxmangl = 3. * gdat.maxmgang
    else:
        if gdat.exprtype == 'sdyn':
            gdat.maxmangl = 1.
        if gdat.exprtype == 'ferm':
            gdat.maxmangl = 15. / gdat.anglfact
        if gdat.exprtype == 'chan':
            gdat.maxmangl = 10. / gdat.anglfact
        if gdat.exprtype == 'hubb':
            gdat.maxmangl = 1. / gdat.anglfact
    retr_axis(gdat, 'anglhalf', 0., gdat.maxmgangdata, gdat.numbanglhalf)
    retr_axis(gdat, 'angl', 0., gdat.maxmangl, gdat.numbangl)
    retr_axis(gdat, 'anglfull', 0., 3. * gdat.maxmgang, gdat.numbanglfull)
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
    # element indices in each population
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

    gdat.listnamechro = ['totl', 'type', 'prop', 'diag', 'save', 'plot', 'proc', 'lpri', 'llik']
    if gdat.elemtype == 'lght':
        gdat.listnamechro += ['lghtpntsprep', 'lghtpntseval', 'lghtmodl']
    gdat.listnamechro += ['expo']
    gdat.listlegdchro = ['Total', 'Type', 'Proposal', 'Diagnostics', 'Save', 'Plot', 'Process', 'Prior', 'Posterior']
    if gdat.elemtype == 'lght':
        gdat.listlegdchro += ['PS Prep', 'PSF Evaluation', 'Emission model']
    gdat.listlegdchro += ['Exposure']
    gdat.indxchro = dict()
    for k, name in enumerate(gdat.listnamechro):
        gdat.indxchro[name] = k
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
        retr_axis(gdat, 'defl', -gdat.maxmgang, gdat.maxmgang, 40)
        retr_axis(gdat, 'deflelem', -gdat.maxmgang * 1e-2, gdat.maxmgang * 1e-2, 40)

    # lensing problem setup
    ## number of deflection components to plot

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
    gdat.meanlgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
    gdat.meanbgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
    gdat.meanlgalcartmesh, gdat.meanbgalcartmesh = meshgrid(gdat.meanlgalcart, gdat.meanbgalcart)
    if gdat.pixltype == 'cart':
        gdat.apix = (2. * gdat.maxmgangdata / gdat.numbsidecart)**2
        gdat.sizepixl = sqrt(gdat.apix)
        gdat.indxpixlrofi = arange(gdat.numbsidecart**2)
        gdat.indxsidecart = arange(gdat.numbsidecart)
        gdat.indxsidemesh = meshgrid(gdat.indxsidecart, gdat.indxsidecart, indexing='ij')
        gdat.lgalgrid = gdat.meanlgalcart[gdat.indxsidemesh[0].flatten()]
        gdat.bgalgrid = gdat.meanbgalcart[gdat.indxsidemesh[1].flatten()]
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
            # determine the maximum angle at which the contribution of the element will be computed
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
        gdat.numbpixlprox = zeros(gdat.numbprox) 
        for k in range(gdat.numbprox):
            for m in range(len(gdat.indxpixlprox[k])):
                gdat.numbpixlprox[k] += len(gdat.indxpixlprox[k][m])
            gdat.numbpixlprox[k] /= len(gdat.indxpixlprox[k])
        print 'gdat.numbpixlprox'
        print gdat.numbpixlprox
    else:
        gdat.numbpixlprox = gdat.numbpixl

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
        if gdat.fittnumbtrap > 0 and k in gdat.fittindxfixpnumbpnts:
            thisbool = False
        else:
            if gdat.fittnumbtrap > 0:
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
    gdat.propfixp = gdat.propmeanpnts or gdat.propdist or gdat.propbacp or gdat.proppsfp or gdat.proplenp
    if not gdat.propfixp and not gdat.propcomp:
        raise Exception('Either a fixed dimensional parameter or an element component must be perturbed.')
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
    if gdat.fittnumbtrap > 0:
        gdat.numbstdp = gdat.numbfixpprop + gdat.fittmaxmnumbcomp
    else:
        gdat.numbstdp = gdat.numbfixpprop
    gdat.numbstdpfixp = gdat.numbfixpprop
    if gdat.indxfixpprop.size > 0:
        gdat.strgstdp = concatenate((array(gdat.fittlablfixp)[gdat.indxfixpprop], gdat.fittliststrgcomptotl))
        gdat.namestdp = concatenate((array(gdat.fittnamefixp)[gdat.indxfixpprop], gdat.fittliststrgcomptotl))
    else:
        gdat.strgstdp = gdat.fittliststrgcomptotl
        gdat.namestdp = gdat.fittliststrgcomptotl
    gdat.strgstdp = list(gdat.strgstdp)
    gdat.indxstdp = arange(gdat.numbstdp)
    gdat.indxstdpfixp = arange(gdat.numbfixpprop)
    gdat.indxstdpcomp = setdiff1d(gdat.indxstdp, gdat.indxstdpfixp)
    gdat.indxparaprop = gdat.indxfixpprop

    # proposal scale indices for each parameter
    indxelemfull = [range(gdat.fittmaxmnumbpnts[l]) for l in gdat.fittindxpopl]
    indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, 'fitt')
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
            indxstdp = gdat.indxstdppara[valu]
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
            if gdat.propdist:
                gdat.stdvstdp[gdat.fittindxfixpdist+gdat.fittnumbpopl] = 1e-2
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
            print 'fixp'
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
            
            print 'elem'
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
    gdat.indxproptypecomp = []
    if gdat.propwithsing:
        for k in gdat.indxstdp:    
            gdat.indxproptypewith = cntr.incr()
            gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{W}_{%d}$' % k)
            gdat.legdproptype = append(gdat.legdproptype, gdat.fittnamepara[where(gdat.indxstdppara == k)[0][0]])
            gdat.nameproptype = append(gdat.nameproptype, gdat.namestdp[k])
            if k >= gdat.indxstdplgal:
                gdat.indxproptypecomp.append(gdat.indxproptypewith)
    else:    
        gdat.indxproptypewith = cntr.incr()
        gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{W}$')
        gdat.legdproptype = append(gdat.legdproptype, 'Within-model')
        gdat.nameproptype = append(gdat.nameproptype, 'with')
    
    if gdat.fittnumbtrap > 0.:
        
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
    
    if gdat.datatype == 'mock':
        for name in gdat.truenamefixp:
            setattr(gdat, 'true' + name, gdat.truesampvarb[getattr(gdat, 'trueindxfixp' + name)])

    # for each parameter in the fitting model, determine if there is a corresponding parameter in the generative model
    for k in gdat.indxvarbscal:
        try:
            temp = getattr(gdat, 'true' + gdat.listnamevarbscal[k])
        except:
            temp = None
        setattr(gdat, 'corr' + gdat.listnamevarbscal[k], temp)

	gdat.fittcorrfixp = empty(gdat.fittnumbfixp)
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
                    pdfnpriotemp = empty((gdat.numbsidecart + 1, gdat.numbsidecart + 1))
                    for k in range(gdat.numbsidecart + 1):
                        for n in range(gdat.numbsidecart + 1):
                            # temp
                            temp = retr_rele(gdat, gdat.truelenscnts[0, :, :, 0], gdat.binslgalcart[k], gdat.binsbgalcart[n], 1., \
                                                                                                    gdat.trueascadistmeanpop0, gdat.trueacutdistmeanpop0, gdat.indxpixl)
    
                            #temp /= amax(temp)
                            pdfnpriotemp[k, n] = temp**gdat.relnpowr
                    lpdfprio, lpdfprioobjt = retr_spatprio(gdat, pdfnpriotemp)
                    lpdfpriointp = lpdfprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
        setattr(gdat, strgmodl + 'lpdfprio', lpdfprio)
        setattr(gdat, strgmodl + 'lpdfprioobjt', lpdfprioobjt)
        setattr(gdat, strgmodl + 'lpdfpriointp', lpdfpriointp)
        
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

    if gdat.datatype == 'inpt':
        for strgfeat in gdat.fittliststrgfeattotl:
            hist = empty((1, gdat.numbbinsplot))
            bins = getattr(gdat, 'bins' + strgfeat)
            # temp -- this takes all reference elements inside the ROI
            comp = getattr(gdat, 'expr' + strgfeat)
            try:
                hist[0, :] = histogram(comp, bins)[0]
            except:
                hist = None
            setattr(gdat, 'truehist' + strgfeat, hist)


def retr_gradmaps(gdat, maps):
    
    grad = dstack((gradient(maps, gdat.sizepixl, axis=0), gradient(maps, gdat.sizepixl, axis=1)))

    return grad


def retr_rele(gdat, maps, lgal, bgal, defs, asca, acut, indxpixleval, absv=True):
    
    grad = retr_gradmaps(gdat, maps)
        
    defl = retr_defl(gdat, lgal, bgal, defs, 0., 0., asca=asca, acut=acut, indxpixleval=indxpixleval).reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    dotstemp = sum(grad * defl, 2)
    if absv:
        dotstemp = fabs(dotstemp)
    else:
        dotstemp = dotstemp
    
    dots = mean(dotstemp) * gdat.sizepixl
    
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
    

def retr_fromgdat(gdat, gdatmodi, strg, strgvarb, mometype='medi'):
    
    if strgvarb == 'datacnts':
        varb = gdat.datacnts
    elif strg == 'this':
        if mometype == 'errr':
            varb = getattr(gdatmodi, strg + 'errr' + strgvarb)
        else:
            varb = getattr(gdatmodi, strg + strgvarb)
    elif strg == 'post':
        varb = getattr(gdat, mometype + strgvarb)
    else:
        varb = getattr(gdat, strg + strgvarb)

    return copy(varb)


def retr_indxsamp(gdat, strgmodl='fitt'):
    
    spectype = getattr(gdat, strgmodl + 'spectype')
    
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
    
    # features which correlate with significance
    liststrgfeatsign = [[] for l in indxpopl]
    for l in indxpopl:
        for namefeat in gdat.listnamefeatsele:
            liststrgfeatsign[l] += [namefeat]
        
    liststrgcomp = [[] for l in indxpopl]
    listscalcomp = [[] for l in indxpopl]
    for l in indxpopl:
        
        liststrgcomp[l] = ['lgal', 'bgal']
        
        if gdat.elemtype == 'lght':
            liststrgcomp[l] += ['flux']
        if gdat.elemtype == 'lens':
            liststrgcomp[l] += ['defs']
        if gdat.elemtype == 'clus':
            liststrgcomp[l] += ['nobj']
        listscalcomp[l] += ['powrslop']
        
        if gdat.numbener > 1 and gdat.elemtype == 'lght':
            liststrgcomp[l] += ['sind']
            listscalcomp[l] += ['gaus']
            if spectype[l] == 'curv':
                liststrgcomp[l] += ['curv']
                listscalcomp[l] += ['gaus']
            if spectype[l] == 'expo':
                liststrgcomp[l] += ['expo']
                listscalcomp[l] += ['gaus']
        if gdat.elemtype == 'lens':
            if gdat.variasca:
                liststrgcomp[l] += ['asca']
            if gdat.variacut:
                liststrgcomp[l] += ['acut']
        
    listscalcomp = [[] for l in indxpopl]
    for l in indxpopl:
        listscalcomp[l] = ['self', 'self', 'powrslop']
        if gdat.elemtype == 'lght' and gdat.numbener > 1:
            listscalcomp[l] += ['self']
            if spectype[l] == 'curv':
                listscalcomp[l] += ['self']
            if spectype[l] == 'expo':
                listscalcomp[l] += ['self']
        if gdat.elemtype == 'lens':
            if gdat.variacut:
                listscalcomp[l] += ['self']
            if gdat.variacut:
                listscalcomp[l] += ['self']
    
    # variables for which whose marginal distribution and pair-correlations will be plotted
    liststrgfeatodim = [[] for l in indxpopl]
    for l in indxpopl:
        liststrgfeatodim[l] = deepcopy(liststrgcomp[l])
        liststrgfeatodim[l] += ['deltllik', 'gang', 'aang']
        if gdat.elemtype == 'lens':
            liststrgfeatodim[l] += ['mcut', 'diss', 'rele', 'reln', 'reld', 'relc']
        # temp
        if strgmodl == 'true':
            for namefeat in gdat.listnamefeatsele:
                for namesele in gdat.listnamesele:
                    liststrgfeatodim[l] += [namefeat + namesele]
 
    # variables for which pair-correlations will be plotted
    liststrgfeatcorr = [[] for l in indxpopl]
    if gdat.plotelemcorr:
        for l in indxpopl:
            for strgfeat in liststrgfeatodim[l]:
                liststrgfeatcorr[l].append(strgfeat)
    
    # defaults
    liststrgpdfnmodu = [[] for l in indxpopl]
    liststrgfeatmodu = [[] for l in indxpopl]
    for l in indxpopl:
        if gdat.elemtype == 'lght': 
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
        if gdat.elemtype == 'lght':
            liststrgfeat[l] += ['cnts', 'spec', 'specplot']
        if gdat.elemtype == 'lens':
            liststrgfeat[l] += ['deflprof']
    
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
    
    numbdeflsubhplot = 2
    numbdeflsingplot = numbdeflsubhplot
    if numbtrap > 0:
        numbdeflsingplot += 3

    cntr = tdpy.util.cntr()
    
    liststrgfeatdefa = deepcopy(liststrgfeattotl)
    if gdat.elemtype == 'lght':
        for strgfeat in ['sind', 'curv', 'expo']:
            if not strgfeat in liststrgfeatdefa:
                liststrgfeatdefa.append(strgfeat)

    dicttemp = {}
    if numbtrap > 0:

        # number of elements
        for l in indxpopl:
            dicttemp['indxfixpnumbpntspop%d' % l] = cntr.incr()
        
        # hyperparameters
        ## mean number of elements
        for l in indxpopl:
            dicttemp['indxfixpmeanpntspop%d' % l] = cntr.incr()
        
        for strgfeat in liststrgfeatpriototl:
            try:
                disttype = getattr(gdat, strgmodl + strgfeat + 'disttype')
            except:
                if strgfeat == 'lgal' or strgfeat == 'bgal':
                    disttype = 'self'
                setattr(gdat, strgmodl + strgfeat + 'disttype', disttype)

        ## distribution shapes
        liststrgvarb = []
        for l in indxpopl:
            for pdfnfeat, strgfeat in zip(liststrgpdfnprio[l], liststrgfeatprio[l]):
                if pdfnfeat == 'exposcal':
                    liststrgvarb += [strgfeat + 'distscal']
                if pdfnfeat == 'powrslop':
                    liststrgvarb += [strgfeat + 'distslop']
                if pdfnfeat == 'gausmean':
                    liststrgvarb += [strgfeat + 'distmean']
                if pdfnfeat == 'gausstdv':
                    liststrgvarb += [strgfeat + 'diststdv']
                if pdfnfeat == 'gausmeanstdv':
                    liststrgvarb += [strgfeat + 'distmean', strgfeat + 'diststdv']
        # temp
        for strgvarb in liststrgvarb:
            strgtemp = 'indxfixp' + strgvarb
            temp = zeros(numbpopl, dtype=int) - 1
            dicttemp[strgtemp] = temp

        for l in indxpopl:
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
                if liststrgpdfnprio[l][k] == 'gausmean' or liststrgpdfnprio[l][k] == 'gausmeanstdv':
                    dicttemp['indxfixp' + strgfeatprio + 'distmeanpop%d' % l] = cntr.incr()
                    dicttemp['indxfixp' + strgfeatprio + 'distmean'][l] = dicttemp['indxfixp' + strgfeatprio + 'distmeanpop%d' % l]
                if liststrgpdfnprio[l][k] == 'gausstdv' or liststrgpdfnprio[l][k] == 'gausmeanstdv':
                    dicttemp['indxfixp' + strgfeatprio + 'diststdvpop%d' % l] = cntr.incr()
                    dicttemp['indxfixp' + strgfeatprio + 'diststdv'][l] = dicttemp['indxfixp' + strgfeatprio + 'diststdvpop%d' % l]
    
        dicttemp['indxfixpnumbpnts'] = []
        dicttemp['indxfixpmeanpnts'] = []
        for strg, valu in dicttemp.iteritems():
            if strg[8:].startswith('numbpntsp'):
                dicttemp['indxfixpnumbpnts'].append(valu)
            if strg[8:].startswith('meanpntsp'):
                dicttemp['indxfixpmeanpnts'].append(valu)
        dicttemp['indxfixpnumbpnts'] = array(dicttemp['indxfixpnumbpnts'])
        dicttemp['indxfixpmeanpnts'] = array(dicttemp['indxfixpmeanpnts'])
    
        dicttemp['indxfixpdist'] = []
        for strg, valu in dicttemp.iteritems():
            if strg[12:16] == 'dist' and isscalar(valu):
                dicttemp['indxfixpdist'].append(valu)
            
        dicttemp['indxfixpdist'] = array(dicttemp['indxfixpdist']) 
        dicttemp['indxfixphypr'] = array(dicttemp['indxfixpdist'] +  dicttemp['indxfixpmeanpnts'])
    
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            if psfntype == 'singgaus':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
            if psfntype == 'doubking':
                dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpsigtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixpgamtene%devt%d' % (i, m)] = cntr.incr()
                dicttemp['indxfixppsffene%devt%d' % (i, m)] = cntr.incr()
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
    scalmeanpnts = getattr(gdat, strgmodl + 'scalmeanpnts')

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
                scalfixp[k] = scalmeanpnts
    
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
                if gdat.elemtype == 'lens' and gdat.fittstdvdefsdistslop != 'none':
                    scalfixp[k] = 'gaus'
                else:
                    scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('sinddistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablsind)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('sinddiststdv'):
                lablfixp[k] = r'$\sigma_{%s%s}$' % (strgpoplcomm, gdat.lablsind)
                scalfixp[k] = 'logt'
            
            if namefixp[k].startswith('curvdistmean'):
                lablfixp[k] = r'$\lambda_{%s%s}$' % (strgpoplcomm, gdat.lablcurv)
                scalfixp[k] = 'logt'
            
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
                    scalfixp[k] = 'logt'
                if strgvarb.startswith('psff'):
                    scalfixp[k] = 'self'
                if strgvarb.startswith('onor'):
                    scalfixp[k] = 'logt'
                if strgvarb.startswith('oind'):
                    scalfixp[k] = 'logt'
            
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


def setp_namevarbsing(gdat, strgvarb, popl, ener, evtt, back):
    
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


def setp_namevarbvalu(gdat, strgvarb, valu, popl=False, ener=False, evtt=False, back=False, strgmodl='true'):
    
    liststrgvarb = setp_namevarbsing(gdat, strgvarb, popl, ener, evtt, back)
    
    for strgvarb in liststrgvarb:
        try:
            strgvarbtemp = getattr(gdat, strgmodl + strgvarb)
            if strgvarbtemp == None:
                raise
        except:
            setattr(gdat, strgmodl + strgvarb, valu)
 

def setp_namevarblimt(gdat, strgvarb, listvalu, typelimt='minmmaxm', popl=False, ener=False, evtt=False, back=False, strgmodl='true'):
    
    liststrgvarb = setp_namevarbsing(gdat, strgvarb, popl, ener, evtt, back)

    for strgvarb in liststrgvarb:

        if typelimt == 'minmmaxm':
            try:
                minmpara = getattr(gdat, strgmodl + 'minm' + strgvarb)
                if minmpara == None:
                    raise
            except:
                setattr(gdat, strgmodl + 'minm' + strgvarb, listvalu[0])
            try:
                maxmpara = getattr(gdat, strgmodl + 'maxm' + strgvarb)
                if maxmpara == None:
                    raise
            except:
                setattr(gdat, strgmodl + 'maxm' + strgvarb, listvalu[1])
        else:
            try:
                meanpara = getattr(gdat, strgmodl + 'mean' + strgvarb)
                if meanpara == None:
                    raise
            except:
                setattr(gdat, strgmodl + 'mean' + strgvarb, listvalu[0])
            try:
                stdvpara = getattr(gdat, strgmodl + 'stdv' + strgvarb)
                if stdvpara == None:
                    raise
            except:
                setattr(gdat, strgmodl + 'stdv' + strgvarb, listvalu[1])
 

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
    
    # temp
    #if strg == 'this' or gdatmodi != None:
    #    pathfold = gdat.pathplot + gdat.namesampdist + '/fram/'
    #elif strg == 'true' or strg == '':
    #    pathfold = gdat.pathinit
    #elif strg == 'post':
    #    pathfold = gdat.pathplot + gdat.namesampdist + '/finl/'
    
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
    
    nameplot = '%s%s%s%s%s' % (strgplot, strgener, strgevtt, strgpopl, strgswep)
   
    path = retr_plotpath(gdat, gdatmodi, strg, nameplot)
    
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
    
    # transdimensional elements
    if strg == 'post' or strg == 'this':
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
    
    # fixed-dimensional objects
    if gdat.elemtype == 'lens':
        if strg == 'this':
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
                mrkrsize = retr_mrkrsize(gdat, getattr(gdat, 'true' + gdat.namecompampl)[l][0, :])
                lgal = copy(gdat.truelgal[l][0, :])
                bgal = copy(gdat.truebgal[l][0, :])
                numbpnts = int(gdat.truenumbpnts[l])
                
                if gdatmodi != None:
                   
                    ## associations
                    ### missed
                    indx = gdatmodi.thistrueindxpntsasscmiss[l]
                    axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphpnts, label=gdat.legdtruemiss, facecolor='none', \
                                                                                                                                marker='s', linewidth=2, color='g')
                    
                    ### hit
                    indx = gdatmodi.thistrueindxpntsasschits[l]
                    
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
                mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][l]])
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
        ampl = zeros(gdat.numbprvlhigh)
        cntr = 0
        for r in gdat.indxstkscond:
            if r in gdat.indxprvlhigh:
                lgal[cntr] = gdat.dictglob['poststkscond'][r]['lgal'][0]
                bgal[cntr] = gdat.dictglob['poststkscond'][r]['bgal'][0]
                ampl[cntr] = gdat.dictglob['poststkscond'][r][gdat.namecompampl][0]
                cntr += 1
        mrkrsize = retr_mrkrsize(gdat, ampl)
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
   
    #fluxpare, lgalpare, bgalpare, sindpare, fluxauxi, radiauxi, anglauxi, sindauxi \
    #                                                        = sympy.symbols('fluxpare lgalpare bgalpare sindpare fluxauxi radiauxi anglauxi sindauxi')
    
    fluxpare, lgalpare, bgalpare, fluxauxi, lgalauxi, bgalauxi = sympy.symbols('fluxpare lgalpare bgalpare fluxauxi lgalauxi bgalauxi')
    
    #matr = sympy.Matrix([[     fluxauxi, 0, 0 , 0,                         fluxpare,                                    0,                                               0, 0], \
    #                     [            0, 1, 0 , 0, -radiauxi * sympy.cos(anglauxi) , (1 - fluxauxi) * sympy.cos(anglauxi), (fluxauxi - 1) * radiauxi * sympy.sin(anglauxi), 0], \
    #                     [            0, 0, 1 , 0, -radiauxi * sympy.cos(anglauxi) , (1 - fluxauxi) * sympy.cos(anglauxi), (1 - fluxauxi) * radiauxi * sympy.cos(anglauxi), 0], \
    #                     [            0, 0, 0 , 1,                                0,                                    0,                                               0, 0], \
    #                     [ 1 - fluxauxi, 0, 0 , 0,                        -fluxpare,                                    0,                                               0, 0], \
    #                     [            0, 1, 0 , 0, -radiauxi * sympy.cos(anglauxi) , -fluxauxi * sympy.cos(anglauxi), fluxauxi * radiauxi * sympy.sin(anglauxi), 0], \
    #                     [            0, 0, 1 , 0, -radiauxi * sympy.cos(anglauxi) , -fluxauxi * sympy.cos(anglauxi), -fluxauxi * radiauxi * sympy.cos(anglauxi), 0], \
    #                     [            0, 0, 0 , 0,                               0 ,                                    0,                                               0, 1]])

    matr = sympy.Matrix([[     fluxauxi,  fluxpare, 0,            0, 0,            0], \
                         [ 1 - fluxauxi, -fluxpare, 0,            0, 0,            0], \
                         [            0, -lgalauxi, 1, 1 - fluxauxi, 0,            0], \
                         [            0, -lgalauxi, 1,    -fluxauxi, 0,            0], \
                         [            0, -bgalauxi, 0,            0, 1, 1 - fluxauxi], \
                         [            0, -bgalauxi, 0,            0, 1,    -fluxauxi]])

    jcbn = matr.det()
    print jcbn

    return jcbn

#retr_jcbn()

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
    
    deflcutf = defs / fracanglasca
    
    # second term in the NFW deflection profile
    fact = ones_like(fracanglasca)
    indxlowr = where(fracanglasca < 1.)[0]
    indxuppr = where(fracanglasca > 1.)[0]
    fact[indxlowr] = arccosh(1. / fracanglasca[indxlowr]) / sqrt(1. - fracanglasca[indxlowr]**2)
    fact[indxuppr] = arccos(1. / fracanglasca[indxuppr]) / sqrt(fracanglasca[indxuppr]**2 - 1.)
    
    if asym:
        deflcutf *= log(fracanglasca / 2.) + fact
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
def retr_defl_jitt(indxpixl, lgalgrid, bgalgrid, lgal, bgal, bein, ellp, angl, rcor=0., asca=None, acut=None, indxpixleval=None):
    
    if indxpixleval == None:
        indxpixleval = indxpixl
    
    # translate the grid
    lgaltran = lgalgrid[indxpixleval] - lgal
    bgaltran = bgalgrid[indxpixleval] - bgal
    
    if acut != None:
        angl = sqrt(lgaltran**2 + bgaltran**2)
        defl = retr_deflcutf(angl, bein, asca, acut)
        
        defllgal = lgaltran / angl * defl
        deflbgal = bgaltran / angl * defl

    else:
        # rotate the grid
        lgalrttr = cos(angl) * lgaltran - sin(angl) * bgaltran
        bgalrttr = sin(angl) * lgaltran + cos(angl) * bgaltran
        
        axisrati = 1. - ellp
        facteccc = sqrt(1. - axisrati**2)
        factrcor = sqrt(axisrati**2 * lgalrttr**2 + bgalrttr**2)
        defllgalrttr = bein * axisrati / facteccc *  arctan(facteccc * lgalrttr / factrcor)
        deflbgalrttr = bein * axisrati / facteccc * arctanh(facteccc * bgalrttr / factrcor)
        
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
    
    # common dictionary
    dicttemp = {}
           
    # temp
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        indxfixpmeanpnts = getattr(gdat, strgmodl + 'indxfixpmeanpnts')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        if gdat.elemtype == 'lght':
            minmflux = getattr(gdat, strgmodl + 'minmflux')
        if gdat.elemtype == 'lens':
            minmdefs = getattr(gdat, strgmodl + 'minmdefs')
        if gdat.elemtype == 'clus':
            minmnobj = getattr(gdat, strgmodl + 'minmnobj')
            binsnobj = getattr(gdat, strgmodl + 'binsnobj')
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
    
    # grab the sample vector
    sampvarb = getattr(gdatobjt, strg + 'sampvarb')
    
    psfp = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfp')]
    bacp = sampvarb[getattr(gdat, strgmodl + 'indxfixpbacp')]
    
    if numbtrap > 0:
        indxelemfull = list(getattr(gdatobjt, strg + 'indxelemfull'))
        # temp -- this may slow down execution
        indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, strgmodl)
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
                        raise Exception('')
                if gdat.variacut:
                    indx = where(sampvarb[indxsampcomp['asca'][l]] < 0.)[0]
                    if indx.size > 0:
                        print 'Asca went negative'
                        sampvarb[indxsampcomp['asca'][l]][indx] = 1e-3 * gdat.anglfact
                        #raise Exception('')

        if gdat.elemtype == 'lght':
            for l in range(numbpopl):
                dicttemp['spec'][l] = retr_spec(gdat, dicttemp['flux'][l], dicttemp['sind'][l], dicttemp['curv'][l], dicttemp['expo'][l], spectype=spectype[l])
        
        for strgfeat in gdat.liststrgfeatconc:
            if strgfeat == 'spec':
                dicttemp['specconc'] = concatenate(dicttemp['spec'], axis=1)
            else:
                dicttemp[strgfeat + 'conc'] = concatenate(dicttemp[strgfeat])
        
        numbpntsconc = dicttemp['lgalconc'].size
    else:
        numbpntsconc = 0
    
    if gdat.verbtype > 1:
        for l in indxpopl:
            for strgcomp in liststrgcomp[l]:
                print strgcomp
                print dicttemp[strgcomp][l]
        print

    # log-prior
    initchro(gdat, gdatmodi, 'lpri')

    lpri = zeros(gdat.numblpri)
    if numbtrap > 0:
        
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
    
        meanpnts = sampvarb[indxfixpmeanpnts]
        
        for l in gdat.fittindxpopl:
            lpri[0] -= 0.5 * gdat.priofactdoff * numbcomp[l] * numbpnts[l]
            lpri[1+0*numbpopl+l] = retr_probpois(numbpnts[l], meanpnts[l])
            
            for k, (strgfeat, pdfnfeat) in enumerate(zip(liststrgfeatprio[l], liststrgpdfnprio[l])):
                
                if pdfnfeat == 'tmpl':

                    if pdfnfeat.endswith('cons'):
                        pdfnspatpriotemp = getattr(gdat, strgmodl + 'pdfnspatpriotemp')
                        spatdistcons = sampvarb[getattr(gdat, strgmodl + 'indxfixpspatdistcons')]
                        lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons)
                        lpdfspatpriointp = lpdfspatprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
                        
                        # temp
                        lpdfspatpriointp = lpdfspatpriointp.T
                        
                        setattr(gdatobjt, strg + 'lpdfspatpriointp', lpdfspatpriointp)
                        setattr(gdatobjt, strg + 'lpdfspatprioobjt', lpdfspatprioobjt)
            
                    else:
                        lpdfspatprioobjt = gdat.fittlpdfspatprioobjt
                
                if pdfnfeat == 'self':
                    minmfeat = getattr(gdat, 'minm' + strgfeat)
                    maxmfeat = getattr(gdat, 'maxm' + strgfeat)
                    # temp -- this may be sped up a bit
                    lpri[1+(k+1)*numbpopl+l] = numbpnts[l] * log(1. / (maxmfeat - minmfeat))
                
                if False:
                    if strgpdfn == 'disc':
                        indxfixpbgaldistscal = getattr(gdat, strgmodl + 'indxfixpbgaldistscalpop%d' % l)
                        lpri[1+2*numbpopl+l] = sum(log(pdfn_dexp(dicttemp['bgal'][l], gdat.maxmgang, sampvarb[indxfixpbgaldistscal]))) 
                    elif strgpdfn == 'exposcal':
                        gang = retr_gang(dicttemp['lgal'][l], dicttemp['bgal'][l])
                        indxfixpgangdistscal = getattr(gdat, strgmodl + 'indxfixpgangdistscalpop%d' % l)
                        lpri[1+numbpopl+l] = sum(log(pdfn_expo(gang, gdat.maxmgang, sampvarb[indxfixpgangdistscal]))) 
                        lpri[1+2*numbpopl+l] = -numbpnts[l] * log(2. * pi) 
                    elif strgpdfn == 'tmpl':
                        lpri[1+numbpopl+l] = sum(lpdfspatprioobjt(dicttemp['bgal'][l], dicttemp['lgal'][l], grid=False))
                
                if pdfnfeat == 'powrslop':
                    lpri[1+(k+1)*numbpopl+l] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, dicttemp[strgfeat][l], strgfeat, sampvarb, l)
                if pdfnfeat == 'gausmeanstdv':
                    lpri[1+(k+1)*numbpopl+l] = retr_lprigausdist(gdat, gdatmodi, strgmodl, dicttemp[strgfeat][l], strgfeat, sampvarb, l)
            
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
                        gdatmodi.thislpau[k] = -log(2. * gdat.fittmaxmgang)
                    if listscalcomp[gdatmodi.indxpoplmodi][k] == 'self':
                        minm = getattr(gdat, 'minm' + strgcomp)
                        maxm = getattr(gdat, 'maxm' + strgcomp)
                        # temp -- this can be sped up
                        #fact = getattr(gdat, 'fact' + strgcomp)
                        gdatmodi.thislpau[k] = -log(maxm - minm)
                    if listscalcomp[gdatmodi.indxpoplmodi][k] == 'powrslop':
                        gdatmodi.thislpau[k] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                           strgcomp, gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                    if listscalcomp[gdatmodi.indxpoplmodi][k] == 'gaus':
                        gdatmodi.thislpau[k] = retr_lprigausdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                           strgcomp, gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                
            if gdatmodi.propsplt or gdatmodi.propmerg:
            
                for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi]):
                    if k == 0 or k == 1:
                        gdatmodi.thislpau[gdat.fittmaxmnumbcomp*l+k] = -log(2. *pi) - log(gdat.radispmr) -0.5 * (gdatmodi.auxipara[k] / gdat.radispmr)
                    elif k == 2:
                        if gdatmodi.auxipara[k] < 0. or gdatmodi.auxipara[k] > 1.:
                            gdatmodi.thislpau[k] = -inf
                        else:
                            gdatmodi.thislpau[k] = 0.
                    elif listscalcomp[gdatmodi.indxpoplmodi][k] == 'powrslop':
                        gdatmodi.thislpau[k] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                           strgcomp, gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                    elif listscalcomp[gdatmodi.indxpoplmodi][k] == 'gaus':
                        gdatmodi.thislpau[k] = retr_lprigausdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                           strgcomp, gdatmodi.thissampvarb, gdatmodi.indxpoplmodi)
                
            if gdatmodi.propbrth or gdatmodi.propsplt:
                gdatmodi.thislpau *= -1.
            
            if gdat.verbtype > 1:
                print 'gdatmodi.thislpau'
                print gdatmodi.thislpau
                print 'lpri'
                print lpri

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
                    psfnintp.append(interp1d_pick(gdat.binsangl, psfn[:, :, :, p], axis=1))
            else:
                psfnintp = interp1d_pick(gdat.binsangl, psfn, axis=1)
            # temp
            # setattr(gdatobjt, strg + 'psfnintp', psfnintp)
        
        if gdat.elemtype == 'lens':

            ## subhalos
            if numbpntsconc > 0 and not raww:
                if gdat.propwithsing and gdat.pertmodleval and strg == 'next':
                    if gdatmodi.propfixp:
                        numbelemeval = 0
                        deflelem = gdatmodi.thisdeflelem
                    else:
                        numbelemeval = gdatmodi.numbelemeval
                        deflelem = copy(gdatmodi.thisdeflelem)
                        for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                            dicttemp[namecomp + 'eval'] = gdatmodi.dicttemp[namecomp + 'eval']
                else:
                    numbelemeval = numbpntsconc
                    deflelem = zeros((gdat.numbpixl, 2))
                    for namecomp in getattr(gdat, strgmodl + 'liststrgcomptotl'):
                        dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']
                
                for k in range(numbelemeval):
                    indxpixleval = retr_indxpixleval(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], dicttemp['defseval'][k], gdat.evalcirc)

                    if gdat.variasca:
                        asca = dicttemp['ascaeval'][k]
                    else:
                        asca = gdat.ascaglob

                    if gdat.variacut:
                        acut = dicttemp['acuteval'][k]
                    else:
                        acut = gdat.acutglob

                    deflelem[indxpixleval, :] += retr_defl_jitt(gdat.indxpixl, gdat.lgalgrid, gdat.bgalgrid, \
                                                 dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], dicttemp['defseval'][k], 0., 0., \
                                                                                                                asca=asca, acut=acut, indxpixleval=indxpixleval)
                defl += deflelem
            
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
            
            setattr(gdatobjt, strg + 'defl', defl)
            
        if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
            
            # prepare the evaluation array
            initchro(gdat, gdatmodi, 'lghtpntsprep')
            if gdat.propwithsing and gdat.pertmodleval and strg == 'next':
                if gdatmodi.propfixp:
                    pntsflux = gdatmodi.thispntsflux
                    numbelemeval = 0
                else:
                    pntsflux = copy(gdatmodi.thispntsflux)
                    numbelemeval = gdatmodi.numbelemeval
                    dicttemp['sindeval'] = []
                    dicttemp['curveval'] = []
                    dicttemp['expoeval'] = []
                    for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                        dicttemp[namecomp + 'eval'] = gdatmodi.dicttemp[namecomp + 'eval']
                    dicttemp['speceval'] = retr_spec(gdat, dicttemp['fluxeval'], dicttemp['sindeval'], dicttemp['curveval'], \
                                                                                                                dicttemp['expoeval'], spectype[gdatmodi.indxpoplmodi])
            else:
                pntsflux = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
                numbelemeval = numbpntsconc
                for namecomp in ['lgal', 'bgal', 'spec']:
                    dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']
            stopchro(gdat, gdatmodi, 'lghtpntsprep')
            
            # evaluate the model for each perturbative PS
            initchro(gdat, gdatmodi, 'lghtpntseval')
            for k in range(numbelemeval):
                indxpixleval = retr_indxpixleval(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], dicttemp['speceval'][gdat.indxenerfluxdist[0], k], gdat.evalcirc)
                pntsflux[:, indxpixleval, :] += retr_pntsflux(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], dicttemp['speceval'][:, k], psfnintp, oaxitype, indxpixleval)
            setattr(gdatobjt, strg + 'pntsflux', pntsflux)
            stopchro(gdat, gdatmodi, 'lghtpntseval')

            ### calculate model flux map
            initchro(gdat, gdatmodi, 'lghtmodl')
            modlflux = retr_mapslght(gdat, bacp, pntsflux, gdat.indxcube)
            stopchro(gdat, gdatmodi, 'lghtmodl')

        initchro(gdat, gdatmodi, 'expo')

        ### count map
        modlcnts = retr_cntsmaps(gdat, modlflux)
        setattr(gdatobjt, strg + 'modlcnts', modlcnts)
        
        stopchro(gdat, gdatmodi, 'expo')

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
        
        if gdat.diagmode:
            if not (modlcnts > 0.).all():
                raise Exception('Model prediction is not positive-definite.')

            if not isfinite(llik).all():
                raise Exception('Likelihood is not finite.')

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
        
        numbdeflsubhplot = getattr(gdat, strgmodl + 'numbdeflsubhplot')
        numbdeflsingplot = getattr(gdat, strgmodl + 'numbdeflsingplot')
        
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
                    pntsfluxfull = retr_pntsfluxtotl(gdat, dicttemp['lgalconc'], dicttemp['bgalconc'], dicttemp['specconc'], psfnintp, gdat.fittoaxitype, evalcirc='full')
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

            masshostbein = array([gdat.massfrombein * beinhost**2])
            setattr(gdatobjt, strg + 'masshostbein', masshostbein)
            if numbtrap > 0:
                ### sort with respect to deflection at scale radius
                if numbpntsconc > 0:
                    indxpntssortbrgt = argsort(dicttemp[gdat.namefeatsort + 'conc'])[::-1]
                    for strgcomp in liststrgcomptotl:
                        dicttemp[strgcomp + 'sort'] = dicttemp[strgcomp + 'conc'][indxpntssortbrgt][:numbpntsconc]

            deflsing = zeros((gdat.numbpixl, 2, numbdeflsingplot))
            numbdeflsing = min(numbdeflsubhplot, numbpntsconc) + 2
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
                    indxpixleval = retr_indxpixleval(gdat, dicttemp['lgalsort'][k-3], dicttemp['bgalsort'][k-3], dicttemp['defssort'][k-3], gdat.evalcirc)
                    
                    if gdat.variasca:
                        asca = dicttemp['ascasort'][k-3]
                    else:
                        asca = gdat.ascaglob

                    if gdat.variacut:
                        acut = dicttemp['acutsort'][k-3]
                    else:
                        acut = gdat.acutglob

                    deflsing[indxpixleval, :, k] = retr_defl(gdat, dicttemp['lgalsort'][k-3], dicttemp['bgalsort'][k-3], dicttemp['defssort'][k-3], 0., 0., \
                                                                                                                    asca=asca, acut=acut, indxpixleval=indxpixleval)
                    
            deflsing = deflsing.reshape((gdat.numbsidecart, gdat.numbsidecart, 2, numbdeflsingplot))

            # convergence
            ## total
            conv = retr_conv(gdat, defl) 
            setattr(gdatobjt, strg + 'conv', conv)
            
            ### power spectrum
            #### two dimensional
            convpsec = retr_psec(gdat, conv)
            setattr(gdatobjt, strg + 'convpsec', convpsec)
            
            #### one dimensional
            convpsecodim = retr_psecodim(gdat, convpsec) 
            setattr(gdatobjt, strg + 'convpsecodim', convpsecodim)
            
            ## subhalos
            if numbtrap > 0:
                deflelem = deflelem.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
                
                ### convergence
                convelem = retr_conv(gdat, deflelem) 
                setattr(gdatobjt, strg + 'convelem', convelem)
                
                ####  power spectrum
                ##### two dimensional
                convpsecelem = retr_psec(gdat, convelem)
                setattr(gdatobjt, strg + 'convpsecelem', convpsecelem)

                ##### one dimensional
                convpsecelemodim = retr_psecodim(gdat, convpsecelem) 
                setattr(gdatobjt, strg + 'convpsecelemodim', convpsecelemodim)
            
            ### magnification
            magn = 1. / retr_invm(gdat, defl) 
            setattr(gdatobjt, strg + 'magn', magn)
            
            histdefl = histogram(defl, bins=gdat.binsdefl)[0]
            setattr(gdatobjt, strg + 'histdefl', histdefl)
            if numbtrap > 0:
                histdeflelem = histogram(deflelem, bins=gdat.binsdeflelem)[0]
                setattr(gdatobjt, strg + 'histdeflelem', histdeflelem)
            
            deflsingmgtd = sqrt(sum(deflsing**2, axis=2))
            setattr(gdatobjt, strg + 'deflsing', deflsing)
            setattr(gdatobjt, strg + 'deflsingmgtd', deflsingmgtd)
        
        liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
    
        if gdat.elemtype == 'lens':
            lensfluxmean = mean(sum(lensflux, 3), 0)
            
        ## element related
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

                        indxpixleval = retr_indxpixleval(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['defs'][l][k], gdat.evalcirc)
                        defltemp[indxpixleval, :] -= retr_defl(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['defs'][l][k], \
                                                                                                                   0., 0., asca=asca, acut=acut, indxpixleval=indxpixleval)
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
                        indxpixleval = retr_indxpixleval(gdat, dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], dicttemp['flux'][l][k], gdat.evalcirc)
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
                deflmgtd = sqrt(sum(defl**2, axis=2))
                gradlenscnts = retr_gradmaps(gdat, lenscnts[0, :, :, 0]) * gdat.sizepixl
                gradlenscntsmgtd = sqrt(sum(gradlenscnts**2, axis=2))
                gradlenscnts *= gdat.sizepixl
                indx = where(fabs(gradlenscnts) > 1. * gdat.sizepixl)
                gradlenscnts[indx] = sign(gradlenscnts[indx]) * 1. * gdat.sizepixl

                setattr(gdatobjt, strg + 'deflmgtd', deflmgtd)
                setattr(gdatobjt, strg + 'gradlenscnts', gradlenscnts)
                setattr(gdatobjt, strg + 'gradlenscntsmgtd', gradlenscntsmgtd)

                #### distance to the source
                for l in range(numbpopl):
                    dicttemp['diss'][l] = retr_angldist(gdat, dicttemp['lgal'][l],  dicttemp['bgal'][l], lgalsour, bgalsour)
                
                for l in indxpopl:
                    dicttemp['deflprof'][l] = empty((gdat.numbangl, numbpnts[l]))
                    dicttemp['mcut'][l] = empty(numbpnts[l])
                    dicttemp['rele'][l] = empty(numbpnts[l])
                    dicttemp['reln'][l] = empty(numbpnts[l])
                    dicttemp['reld'][l] = empty(numbpnts[l])
                    dicttemp['relc'][l] = empty(numbpnts[l])
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
                        dicttemp['deflprof'][l][:, k] = retr_deflcutf(gdat.meanangl, dicttemp['defs'][l][k], asca, acut)
             
                        ### truncated mass 
                        dicttemp['mcut'][l][k] = retr_mcut(gdat, dicttemp['defs'][l][k], asca, acut)
                        
                        #### dot product with the source flux gradient
                        # temp -- weigh the energy and PSF bins
                        dicttemp['rele'][l][k] = retr_rele(gdat, lenscnts[0, :, :, 0], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], 1., asca, acut, gdat.indxpixl)
                        dicttemp['reln'][l][k] = dicttemp['rele'][l][k] / dicttemp['defs'][l][k] * gdat.sizepixl
                        dicttemp['reld'][l][k] = retr_rele(gdat, gdat.datacntscart[0, :, :, 0], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], 1., asca, acut, gdat.indxpixl)
                        dicttemp['relc'][l][k] = -retr_rele(gdat, lenscnts[0, :, :, 0], dicttemp['lgal'][l][k], dicttemp['bgal'][l][k], 1., asca, acut, gdat.indxpixl, absv=False)
                   
            if strgmodl == 'true':
                for namesele in gdat.listnamesele:
                    indx = [[] for l in indxpopl]
                    for l in indxpopl:
                        if namesele == 'pars':
                            indx[l] = where(dicttemp['deltllik'][l] > 0.5 * numbcomp[l])[0]
                        if namesele == 'nrel':
                            indx[l] = where(dicttemp['reln'][l] > 0.3)[0]
                        for strgfeat in gdat.trueliststrgfeat[l]:
                            if strgfeat in gdat.listnamefeatsele:
                                dicttemp[strgfeat + namesele][l] = zeros(numbpnts[l])
                                dicttemp[strgfeat + namesele][l][indx[l]] = dicttemp[strgfeat][l][indx[l]]
                    setattr(gdat, 'trueindxelem' + namesele, indx)
                
            # temp
            if strgmodl == 'true' and gdat.verbtype > 0:
                for strgfeat in liststrgfeattotl:
                    minm = getattr(gdat, 'minm' + strgfeat)
                    maxm = getattr(gdat, 'maxm' + strgfeat)
                    if strgfeat == 'gang':
                        minm *= gdat.anglfact
                        maxm *= gdat.anglfact
                        dicttemp[strgfeat][l] *= gdat.anglfact
                    for l in indxpopl:
                        if where(minm > dicttemp[strgfeat][l])[0].size > 0 or where(maxm < dicttemp[strgfeat][l])[0].size > 0:
                            print 'Warning: element feature outside the plot limits.'
                            print 'Feature: '
                            print strgfeat
                            print 'Plot minmimum'
                            print minm
                            print 'Plot maxmimum'
                            print maxm
                            print 'Feature minimum'
                            print amin(dicttemp[strgfeat][l])
                            print 'Feature maximum'
                            print amax(dicttemp[strgfeat][l])
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
                        if strgfeat[:-4] in gdat.listnamefeatsele and strgmodl == 'true':
                            indx = intersect1d(getattr(gdat, strgmodl + 'indxelem' + strgfeat[-4:]), indxmodlpntscomp[l])
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
                        strgpdfn = liststrgpdfnprio[l][liststrgfeatprio[l].index(strgfeat)]
                        minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                        maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                        bins = getattr(gdat, 'bins' + strgfeat)
                        delt = getattr(gdat, 'delt' + strgfeat)
                        
                        xdat = getattr(gdat, strgmodl + 'mean' + strgfeat + 'prio')
                        deltprio = getattr(gdat, strgmodl + 'delt' + strgfeat + 'prio')
                        
                        booltemp = False
                        if strgpdfn.startswith('expo') or strgpdfn.startswith('dexp'):
                            scal = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distscal')[l]]
                            if strgpdfn.startswith('expo'):
                                pdfn = pdfn_expo(xdat, maxm, scal)
                            if strgpdfn.startswith('dexp'):
                                pdfn = pdfn_dexp(xdat, maxm, scal)
                            booltemp = True
                        if strgpdfn.startswith('self'):
                            minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
                            maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
                            pdfn = 1. / (maxm - minm) + zeros_like(xdat)
                            booltemp = True
                        # temp 
                        if strgpdfn.startswith('powrslop'):
                            slop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
                            pdfn = pdfn_powr(xdat, minm, maxm, slop)
                            booltemp = True
                        if strgpdfn.startswith('igam'):
                            cutf = getattr(gdat, 'cutf' + strgfeat)
                            pdfn = pdfn_igam(xdat, slop, cutf)
                            booltemp = True
                        if strgpdfn.startswith('gaus'):
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
        
            # more derived parameters
            if gdat.elemtype == 'lens':
                if numbtrap > 0:
                    if strgmodl != 'true': 
                        if gdat.priofactdoff != 0. or (gdat.fittnamefixp[gdat.fittindxfixpmeanpnts] == 'logt').any():
                            # temp
                            indx = where(sum(gdat.truehistmcut, axis=0) > 0)[0]
                            histmcutcorr = empty((numbpopl, gdat.numbbinsplot))
                            for l in indxpopl:
                                histmcutcorr[l, indx] = gdat.truehistmcutpars[l, indx] * dicttemp['histmcut'][l, indx] / gdat.truehistmcut[l, indx]
                            fracsubhcorr = ones(gdat.numbbinsplot)
                            fracsubhcorr[indx] = array([sum(gdat.truehistmcutpars[:, indx] * dicttemp['histmcut'][:, indx] / \
                                                                                            gdat.truehistmcut[:, indx] * gdat.meanmcut[indx]) / masshostbein])
                            setattr(gdatobjt, strg + 'histmcutcorr', histmcutcorr)
                            setattr(gdatobjt, strg + 'fracsubhcorr', fracsubhcorr)
                
                ## total truncated mass of the subhalo as a cross check
                if gdat.variasca:
                    asca = dicttemp['asca'][0]
                else:
                    asca = gdat.ascaglob
                if gdat.variacut:
                    acut = dicttemp['acut'][0]
                else:
                    acut = gdat.acutglob
                factmcutfromdefs = retr_factmcutfromdefs(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour, asca, acut) 
                masssubh = array([sum(factmcutfromdefs * dicttemp['defs'][0])])

                ## calculate the host mass, subhalo mass and its fraction as a function of host halo-centric angle and interpolate at the Einstein radius of the host
                angl = sqrt((gdat.meanlgalcartmesh - lgalhost)**2 + (gdat.meanlgalcartmesh - bgalhost)**2)
                
                if gdat.elemtype == 'lens':
                    listnamevarbmass = []
                    listnamevarbmassscal = ['masshosttotl']
                    listnametemp = ['delt', 'intg']
                    listnamevarbmassvect = []
                    for strgtemp in listnametemp:
                        listnamevarbmassvect.append('masshost' + strgtemp)
                        listnamevarbmassscal.append('masshost' + strgtemp + 'bein')
                    if numbtrap > 0:
                        listnamevarbmassscal.append('masssubhtotl')
                        listnamevarbmassscal.append('fracsubhtotl')
                        for strgtemp in listnametemp:
                            listnamevarbmassvect.append('masssubh' + strgtemp)
                            listnamevarbmassvect.append('fracsubh' + strgtemp)
                            listnamevarbmassscal.append('masssubh' + strgtemp + 'bein')
                            listnamevarbmassscal.append('fracsubh' + strgtemp + 'bein')

                # find the host, subhalo masses and subhalo mass fraction as a function of halo-centric radius
                for name in listnamevarbmassvect:
                    dicttemp[name] = zeros(gdat.numbanglhalf)
                    for k in gdat.indxanglhalf:
                        if name[4:8] == 'host':
                            convtemp = conv
                        if name[4:8] == 'subh':
                            convtemp = convelem
                        
                        if name.endswith('delt'):
                            indxpixl = where((gdat.binsanglhalf[k] < angl) & (angl < gdat.binsanglhalf[k+1]))
                            dicttemp[name][k] = sum(convtemp[indxpixl]) * gdat.critmden * gdat.apix * gdat.adishost**2# / 2. / pi * gdat.deltanglhalf[k]
                        if name.endswith('intg'):
                            indxpixl = where(angl < gdat.meananglhalf[k])
                            dicttemp[name][k] = sum(convtemp[indxpixl]) * gdat.critmden * gdat.apix * gdat.adishost**2# * indxpixl[0].size
                        
                        if name[:4] == 'frac':
                            if dicttemp['masshost' + name[8:]][k] != 0.:
                                dicttemp['fracsubh' + name[8:]][k] = dicttemp['masssubh' + name[8:]][k] / dicttemp['masshost' + name[8:]][k]
                    setattr(gdatobjt, strg + name, dicttemp[name])
                
                # interpolate the host, subhalo masses and subhalo mass fraction at the Einstein radius and save it as a scalar variable
                for name in listnamevarbmassvect:
                    dicttemp[name + 'bein'] = interp(beinhost, gdat.meananglhalf, dicttemp[name])
                    setattr(gdatobjt, strg + name + 'bein', array([dicttemp[name + 'bein']]))
            
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
        if strg != 'true' and gdat.trueinfo:
            
            if gdat.elemtype == 'lens' and gdat.datatype == 'mock':
                resideflsing = deflsing - gdat.truedeflsing
                residefl = defl - gdat.truedefl
                resideflmgtd = sqrt(sum(residefl**2, axis=2))
                resideflperc = 100. * resideflmgtd / gdat.truedeflmgtd
                resideflsingmgtd = sqrt(sum(resideflsing**2, axis=2))
                resideflsingperc = 100. * resideflsingmgtd / gdat.truedeflsingmgtd
                setattr(gdatobjt, strg + 'resideflsing', resideflsing)
                setattr(gdatobjt, strg + 'residefl', residefl)
                setattr(gdatobjt, strg + 'resideflmgtd', resideflmgtd)
                if numbtrap > 0:
                    resiconvelem = convelem - gdat.trueconvelem
                    resiconvelemperc = 100. * resiconvelem / gdat.trueconvelem
                    setattr(gdatobjt, strg + 'resiconvelem', resiconvelem)
                    setattr(gdatobjt, strg + 'resiconvelemperc', resiconvelemperc)
                resimagn = magn - gdat.truemagn
                resimagnperc = 100. * resimagn / gdat.truemagn
                setattr(gdatobjt, strg + 'resimagn', resimagn)
                setattr(gdatobjt, strg + 'resimagnperc', resimagnperc)
    
        if numbtrap > 0:
            
            # correlate the catalog sample with the reference catalog
            # temp
            if strg != 'true' and gdat.trueinfo:
                trueindxpntsasscmiss = [[] for l in gdat.trueindxpopl]
                trueindxpntsasschits = [[] for l in gdat.trueindxpopl]
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
                    
                    for k in range(numbpnts[l]):
                       
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
                            trueindxpntsasscmiss[l].append(k)
                        else:
                            trueindxpntsasschits[l].append(k)
                    temp = where(indxmodlpnts >= 0)[0]
                    featassc[strgfeat][l][temp] = dicttemp[strgfeat][l][indxmodlpnts[temp]]

                    indxfittpntsassc[l] = unique(indxmodlpnts[where(indxmodlpnts > -1)])
                    indxfittpntsfals[l] = setdiff1d(arange(numbpnts[l]), indxfittpntsassc[l])
                    
                    cmpl[l] = array([float(len(trueindxpntsasschits[l])) / gdat.truenumbpnts[l]])
                    fdis[l] = array([float(indxfittpntsfals[l].size) / numbpnts[l]])
                    
                setattr(gdatobjt, strg + 'trueindxpntsasscmiss', trueindxpntsasscmiss)
                setattr(gdatobjt, strg + 'trueindxpntsasschits', trueindxpntsasschits)
                
                for l in gdat.trueindxpopl:
                    setattr(gdatobjt, strg + 'cmplpop%d' % l, cmpl[l])
                
                for l in gdat.fittindxpopl:
                    setattr(gdatobjt, strg + 'fdispop%d' % l, fdis[l])
                
                for strgfeat in liststrgfeatodimtotl:
                    setattr(gdatobjt, strg + strgfeat + 'assc', featassc[strgfeat])
               
                # completeness
                for strgfeat in liststrgfeatodimtotl:
                    try:
                        truehistfeat = getattr(gdat, 'truehist' + strgfeat)
                        cmplfeat = empty((gdat.truenumbpopl, gdat.numbbinsplot))
                        errrcmplfeat = empty((2, gdat.truenumbpopl, gdat.numbbinsplot))
                        truefeat = getattr(gdat, 'true' + strgfeat)
                        for l in gdat.trueindxpopl:
                            indx = where(isfinite(featassc[strgfeat][l]))[0]
                            bins = getattr(gdat, 'bins' + strgfeat)
                            histfeatassc = histogram(truefeat[l][0, trueindxpntsasschits], bins=bins)[0]
                            if truehistfeat != None:
                                cmplfeat[l, :] = histfeatassc / truehistfeat[l, :]
                                errrcmplfeat[:, l, :] = (cmplfeat[l, :] / sqrt(maximum(ones(gdat.numbbinsplot), truehistfeat[l, :])))[None, :]
                        setattr(gdatobjt, strg + 'cmpl' + strgfeat, cmplfeat)
                        setattr(gdatobjt, strg + 'errrcmpl' + strgfeat, errrcmplfeat)
                    except:
                        pass
                    
                # false discovery rate
                for strgfeat in liststrgfeatodimtotl:
                    fdisfeat = empty((gdat.fittnumbpopl, gdat.numbbinsplot))
                    errrfdisfeat = empty((2, gdat.fittnumbpopl, gdat.numbbinsplot))
                    for l in gdat.fittindxpopl:
                        indx = where(isfinite(featassc[strgfeat][l]))[0]
                        bins = getattr(gdat, 'bins' + strgfeat)
                        histfeatfals = histogram(dicttemp[strgfeat][l][indxfittpntsfals[l]], bins=bins)[0]
                        fitthistfeat = getattr(gdatobjt, strg + 'hist' + strgfeat)
                        fdisfeat[l, :] = histfeatfals / fitthistfeat[l, :]
                        errrfdisfeat[:, l, :] = (fdisfeat[l, :] / sqrt(maximum(ones(gdat.numbbinsplot), fitthistfeat[l, :])))[None, :]
                    setattr(gdatobjt, strg + 'fdis' + strgfeat, fdisfeat)
                    setattr(gdatobjt, strg + 'errrfdis' + strgfeat, errrfdisfeat)


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


def retr_sersprof(spec, angl, sizehalf, indx=4):
   
    # this approximation works for 0.5  < indx < 10
    factsers = 1.9992 * indx - 0.3271

    # surface brightness at the half-light radius
    # temp -- this only works for indx == 4!
    sbrthalf = spec / 7.2 / pi / sizehalf**2

    # surface brightness profile
    sersprof = sbrthalf * exp(-factsers * ((angl / sizehalf)**(1. / indx) - 1.))
     
    return sersprof



