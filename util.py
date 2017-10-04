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


## custom random variables, pdfs, cdfs and icdfs

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


def pdfn_self(xdat, minm, maxm):
    
    pdfn = 1. / (maxm - minm)
    
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


def cdfn_gaus(para, meanpara, stdvpara):
   
    cdfn = 0.5  * (1. + sp.special.erf((para - meanpara) / sqrt(2) / stdvpara))
    
    return cdfn


def icdf_gaus(cdfn, meanpara, stdvpara):
    
    para = meanpara + stdvpara * sqrt(2) * sp.special.erfinv(2. * cdfn - 1.)

    return para


def pdfn_gaus(xdat, mean, stdv):
    
    pdfn = 1. / sqrt(2. * pi) / stdv * exp(-0.5 * ((xdat - mean) / stdv)**2)

    return pdfn


def cdfn_lgau(para, mean, stdv):
    
    cdfn = cdfn_gaus(log(para), log(mean), stdv)

    return cdfn


def icdf_lgau(cdfn, mean, stdv):
    
    icdf = exp(icdf_gaus(cdfn, log(mean), stdv))

    return icdf


def pdfn_lgau(xdat, mean, stdv):
    
    pdfn = pdfn_gaus(log(xdat), log(mean), stdv)

    return pdfn


def cdfn_eerr(para, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    tranpara = (para - meanpara) / stdvpara
    cdfnnormpara = 0.5 * (sp.special.erf(tranpara / sqrt(2.)) + 1.)
    cdfn = (cdfnnormpara - cdfnnormminm) / cdfnnormdiff

    return cdfn


def icdf_eerr(cdfn, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    cdfnnormpara = cdfn * cdfnnormdiff + cdfnnormminm
    tranpara = sp.special.erfinv(2. * cdfnnormpara - 1.) * sqrt(2)
    para = tranpara * stdvpara + meanpara
   
    return para


def cdfn_logt(para, minmpara, factpara):

    cdfn = log(para / minmpara) / factpara

    return cdfn


def icdf_logt(cdfn, minmpara, factpara):
    
    para = exp(cdfn * factpara) * minmpara

    return para


def cdfn_atan(para, minmpara, maxmpara):
    
    cdfn = (arctan(para) - arctan(minmpara)) / (arctan(maxmpara) - arctan(minmpara))
    
    return cdfn


def icdf_atan(cdfn, minmpara, maxmpara):

    para = tan((arctan(maxmpara) - arctan(minmpara)) * cdfn + arctan(minmpara))
    
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
    
    if gdat.diagmode:
        if fixpunit == 0:
            print 'Warning. CDF is zero.'

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


def retr_lprbpois(data, modl):
    
    lprb = data * log(modl) - modl - sp.special.gammaln(data + 1)
    
    return lprb
    
        
def retr_sampvarb(gdat, strgmodl, samp, indxsampcomp=None):
    
    sampvarb = zeros_like(samp)
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    
    for k in indxfixp:
        sampvarb[k] = icdf_fixp(gdat, strgmodl, samp[k], k)
        if gdat.diagmode:
            if not isfinite(sampvarb[k]):
                namefixp = getattr(gdat, strgmodl + 'namefixp')
                print 'namefixp[k]'
                print namefixp[k]
                print 'samp[k]'
                print samp[k]
                raise Exception('ICDF is infinite!')
    if indxsampcomp != None:
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
        listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
        retr_sampvarbtrap(gdat, strgmodl, indxsampcomp, gdat.indxregi, indxpopl, liststrgcomp, listscalcomp, samp, sampvarb)

    return sampvarb
    

def retr_sampvarbtrap(gdat, strgmodl, indxsampcomp, indxregi, indxpopl, liststrgcomp, listscalcomp, samp, sampvarb, indxelemfull=slice(None)):
    
    for l in indxpopl:
        for d in indxregi:
            for g in range(len(listscalcomp[l])):
                sampvarb[indxsampcomp[liststrgcomp[l][g]][l][d][indxelemfull]] = icdf_trap(gdat, strgmodl, samp[indxsampcomp[liststrgcomp[l][g]][l][d][indxelemfull]], \
                                                                                                                         sampvarb, listscalcomp[l][g], liststrgcomp[l][g], l, d)


def icdf_trap(gdat, strgmodl, cdfn, sampvarb, scalcomp, strgcomp, l, d):
    
    if scalcomp == 'self' or scalcomp == 'dexp' or scalcomp == 'expo' or scalcomp == 'powrslop':
        minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
        if scalcomp == 'powrslop':
            maxm = getattr(gdat, strgmodl + 'maxm' + strgcomp)
            distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]
            icdf = icdf_powr(cdfn, minm, maxm, distslop)
        else:
            fact = getattr(gdat, strgmodl + 'fact' + strgcomp)
            icdf = icdf_self(cdfn, minm, fact)
    if scalcomp == 'igam':
        distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distslop')[l]]
        cutf = getattr(gdat, 'cutf' + strgcomp)
        icdf = icdf_igam(cdfn, distslop, cutf)
    if scalcomp == 'gaus':
        distmean = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distmean')[l]]
        diststdv = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'diststdv')[l]]
        icdf = icdf_gaus(cdfn, distmean, diststdv)
    
    return icdf


def cdfn_trap(gdat, gdatmodi, icdf):
    
    listscalcomp = gdat.fittlistscalcomp[gdatmodi.indxpoplmodi[0]]
    cdfn = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi[0]])
    for k, strgcomp in enumerate(gdat.fittliststrgcomp[gdatmodi.indxpoplmodi[0]]):
        
        if listscalcomp[k] == 'self' or listscalcomp[k] == 'dexp' or listscalcomp[k] == 'expo' or listscalcomp[k] == 'powrslop':
            minm = getattr(gdat, 'fittminm' + strgcomp)
            if listscalcomp[k] == 'powrslop':
                maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                distslop = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[gdatmodi.indxpoplmodi[0]]]
                cdfn[k] = cdfn_powr(icdf[k], minm, maxm, distslop)
            else:
                fact = getattr(gdat, 'fittfact' + strgcomp)
                cdfn[k] = cdfn_self(icdf[k], minm, fact)
        if listscalcomp[k] == 'igam':
            distslop = gdatmodi.sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[gdatmodi.indxpoplmodi[0]]]
            cutf = getattr(gdat, 'cutf' + strgcomp)
            cdfn[k] = cdfn_igam(icdf[k], distslop, cutf)
        if listscalcomp[k] == 'gaus':
            distmean = gdatmodi.sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[gdatmodi.indxpoplmodi[0]]]
            diststdv = gdatmodi.sampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[gdatmodi.indxpoplmodi[0]]]
            cdfn[k] = cdfn_gaus(icdf[k], distmean, diststdv)
    
    return cdfn


## proposals
### decode the transdimensional element list
def retr_indxsampcomp(gdat, indxelemfull, strgmodl):
    
    numbtrapregipoplcuml = getattr(gdat, strgmodl + 'numbtrapregipoplcuml')
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    boolelemspec = getattr(gdat, strgmodl + 'boolelemspec')
    indxregipopl = getattr(gdat, strgmodl + 'indxregipopl')
    indxcomp = getattr(gdat, strgmodl + 'indxcomp')
    spectype = getattr(gdat, strgmodl + 'spectype')
    indxsampcompinit = getattr(gdat, strgmodl + 'indxsampcompinit')
    liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
    
    indxsampcomp = dict()
    for strgcomp in liststrgcomptotl:
        indxsampcomp[strgcomp] = [[[] for d in indxregipopl[l]] for l in indxpopl]

    indxsampcomp['comp'] = [[[] for d in indxregipopl[l]] for l in indxpopl]
    for l in indxpopl:
        for d in indxregipopl[l]:
            indxsamptemp = indxsampcompinit + numbtrapregipoplcuml[l][d] + array(indxelemfull[l][d], dtype=int) * numbcomp[l]
            if elemtype[l] == 'lghtline':
                indxsampcomp['elin'][l][d] = indxsamptemp
                indxsampcomp['flux'][l][d] = indxsamptemp + 1
            else:
                indxsampcomp['lgal'][l][d] = indxsamptemp
                indxsampcomp['bgal'][l][d] = indxsamptemp + 1
                if boolelemspec[l]:
                    indxsampcomp['flux'][l][d] = indxsamptemp + 2
                    if gdat.numbener > 1:
                        if spectype[l] == 'colr':
                            for i in gdat.indxener:
                                if i == 0:
                                    continue
                                indxsampcomp['sind%04d' % i][l][d] = indxsamptemp + 3 + i - 1 
                        else:
                            indxsampcomp['sind'][l][d] = indxsamptemp + 3
                            if spectype[l] == 'curv':
                                indxsampcomp['curv'][l][d] = indxsamptemp + 4
                            if spectype[l] == 'expc':
                                indxsampcomp['expc'][l][d] = indxsamptemp + 4
                if elemtype[l] == 'lens':
                    indxsampcomp['defs'][l][d] = indxsamptemp + 2
                    if gdat.variasca:
                        indxsampcomp['asca'][l][d] = indxsamptemp + 3
                    if gdat.variacut:
                        indxsampcomp['acut'][l][d] = indxsamptemp + 4
                if elemtype[l] == 'clus':
                    indxsampcomp['nobj'][l][d] = indxsamptemp + 2
            indxsampcomp['comp'][l][d] = repeat(indxsamptemp, numbcomp[l]) + tile(indxcomp[l], len(indxelemfull[l][d]))
    
    if gdat.diagmode:
        numbpara = getattr(gdat, strgmodl + 'numbpara')
        for l in indxpopl:
            for d in indxregipopl[l]:
                if len(indxsampcomp['comp'][l][d]) > 0:
                    if amax(indxsampcomp['comp'][l][d]) > numbpara:
                        print 'ld'
                        print l, d
                        print 'indxsampcomp[comp][l][d]'
                        print indxsampcomp['comp'][l][d]
                        print 'numbpara'
                        print numbpara
                        #print 'sampvarb'
                        #for k in range(sampvarb.size):
                        #    print sampvarb[k]
                        raise Exception('Sample indices of the elements are bad.') 
    
    return indxsampcomp


### update the current state with a Gaussian draw with potential tail
def retr_gaus(gdat, gdatmodi, indxparamodi, stdvpara):
    
    numbparamodi = indxparamodi.size
    stdvtemp = normal(size=numbparamodi) * stdvpara
    if rand() < gdat.fracproprand:
        stdvtemp[choice(indxparamodi)] += randn()
    
    gdatmodi.nextsamp[indxparamodi] = gdatmodi.thissamp[indxparamodi] + stdvtemp
    gdatmodi.nextsamp[indxparamodi] = gdatmodi.nextsamp[indxparamodi] % 1.

       
### update sampler state
def updt_samp(gdat, gdatmodi):
   
    if gdat.verbtype > 1:
        print 'updt_samp()'

    if False:
        print 'gdatmodi.indxsampmodi'
        print gdatmodi.indxsampmodi
        print 'gdatmodi.thissampvarb'
        for k in range(20):
            print gdatmodi.thissampvarb[k]
        print 'gdatmodi.nextsampvarb'
        for k in range(20):
            print gdatmodi.nextsampvarb[k]
    
    # update the sample and the unit sample vectors
    gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
    gdatmodi.thissamp[gdatmodi.indxsampmodi] = gdatmodi.nextsamp[gdatmodi.indxsampmodi]
   
    # update the log-prior
    gdatmodi.thislpritotl = copy(gdatmodi.nextlpritotl)
        
    # update the log-likelihood
    if gdatmodi.propllik:
        for dd, d in enumerate(gdatmodi.indxregieval):
            gdatmodi.thisllik[d][gdatmodi.indxcubeeval[0][dd]] = copy(gdatmodi.nextllik[dd])
        gdatmodi.thislliktotl = copy(gdatmodi.nextlliktotl)

    if gdat.penalpridiff and not gdatmodi.prophypr:
        gdatmodi.thislpridiff = gdatmodi.nextlpridiff
        
    if gdatmodi.proppsfnconv:
        gdatmodi.thispsfnconv = [[[] for i in gdat.indxener] for m in gdat.indxevtt]
        for m in gdat.indxevtt:
            for i in gdat.indxener:
                gdatmodi.thispsfnconv[m][i] = AiryDisk2DKernel(gdatmodi.nextpsfp[i] / gdat.sizepixl)
    
    if gdat.fitthostemistype != 'none':
        if gdatmodi.prophost:
            for dd, d in enumerate(gdatmodi.indxregieval):
                gdatmodi.thissbrthost[d][gdatmodi.indxcubeeval[1][dd]] = copy(gdatmodi.nextsbrthost[dd])
    
    if gdat.fittlensmodltype != 'none':
        if gdatmodi.prophost:
            gdatmodi.thisdeflhost = copy(gdatmodi.nextdeflhost)
    
    if gdatmodi.propelem:
        
        if gdatmodi.proptran:
            gdatmodi.thisindxelemfull = deepcopy(gdatmodi.nextindxelemfull)
            gdatmodi.thisindxsampcomp = deepcopy(gdatmodi.nextindxsampcomp)
        
        # temp
        # keep the state vectors clean
        if gdatmodi.propdeth:
            gdatmodi.thissampvarb[gdatmodi.indxsamptran[0]] = 0.
            gdatmodi.thissamp[gdatmodi.indxsamptran[0]] = 0.
    
        if gdat.fittboolelemfore[gdatmodi.indxpoplmodi[0]]:
            if gdat.verbtype > 1:
                print 'gdatmodi.nextsbrtpnts'
                summgene(gdatmodi.nextsbrtpnts)
                print 'gdatmodi.thissbrtpnts'
                summgene(gdatmodi.thissbrtpnts)
            for dd, d in enumerate(gdatmodi.indxregieval):
                gdatmodi.thissbrtpnts[d][gdatmodi.indxcubeeval[0][dd]] = copy(gdatmodi.nextsbrtpnts[dd])
        if gdat.fittelemtype[gdatmodi.indxpoplmodi[0]] == 'lens':
            for dd, d in enumerate(gdatmodi.indxregieval):
                gdatmodi.thisdeflelem[d] = copy(gdatmodi.nextdeflelem[dd])


def initcompfromstat(gdat, gdatmodi, namerefr):
    
    for l in gdat.fittindxpopl:
        for d in gdat.fittindxregipopl[l]:
            for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
                minm = getattr(gdat, 'fittminm' + strgcomp)
                maxm = getattr(gdat, 'fittmaxm' + strgcomp)
                comp = getattr(gdat, namerefr + strgcomp)[l][d]
                try:
                    if gdat.fittlistscalcomp[l][k] == 'self':
                        fact = getattr(gdat, 'fittfact' + strgcomp)
                        compunit = cdfn_self(comp, minm, fact)
                    if gdat.fittlistscalcomp[l][k] == 'powrslop' or gdat.fittlistscalcomp[l][k] == 'igam':
                        slop = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distslop')[l]]
                        if gdat.fittlistscalcomp[l][k] == 'powrslop':
                            compunit = cdfn_powr(comp, minm, maxm, slop)
                        if gdat.fittlistscalcomp[l][k] == 'igam':
                            cutf = getattr(gdat, 'cutf' + strgcomp)
                            compunit = cdfn_igam(comp, slop, cutf)
                    if gdat.fittlistscalcomp[l][k] == 'gaus':
                        distmean = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'distmean')[l]]
                        diststdv = gdatmodi.thissampvarb[getattr(gdat, 'fittindxfixp' + strgcomp + 'diststdv')[l]]
                        compunit = cdfn_gaus(comp, distmean, diststdv)
                except:
                    if gdat.verbtype > 0:
                        print 'Initialization from the reference catalog failed for %s. Sampling randomly...' % strgcomp
                    compunit = rand(gdat.truenumbelem[l][d])
                gdatmodi.thissamp[gdatmodi.thisindxsampcomp[strgcomp][l][d]] = compunit


## model evaluation

### find the spectra of sources
def retr_spec(gdat, flux, sind=None, curv=None, expc=None, sind0001=None, sind0002=None, elin=None, edis=None, spectype='powr', plot=False):
    
    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if plot:
            meanener = gdat.meanenerplot
        else:
            meanener = gdat.meanener

        if spectype == 'gaus':
            spec = 1. / edis / sqrt(2. * pi) * flux[None, :] * exp(-0.5 * ((gdat.meanener[:, None] - elin[None, :]) / edis)**2)
        if spectype == 'powr':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :])
        if spectype == 'colr':
            if plot:
                spec = zeros((gdat.numbenerplot, flux.size))
            else:
                spec = empty((gdat.numbener, flux.size))
                for i in gdat.indxener:
                    if i == 0:
                        spec[i, :] =  flux
                    if i == 1:
                        spec[i, :] = flux * (gdat.meanener[i] / gdat.enerpivt)**(-sind0001)
                    if i == 2:
                        spec[i, :] = flux * (gdat.meanener[i] / gdat.enerpivt)**(-sind0002)
        if spectype == 'curv':
            spec = flux[None, :] * meanener[:, None]**(-sind[None, :] - gdat.factlogtenerpivt[:, None] * curv[None, :])
        if spectype == 'expc':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :]) * exp(-(meanener - gdat.enerpivt)[:, None] / expc[None, :])
    
    return spec


### find the surface brightness due to one point source
def retr_sbrtpnts(gdat, lgal, bgal, spec, psfnintp, oaxitype, indxpixleval):
    
    # calculate the distance to all pixels from each point source
    dist = retr_angldistunit(gdat, lgal, bgal, indxpixleval)
   
    # interpolate the PSF onto the pixels
    if gdat.kernevaltype == 'ulip':
        if oaxitype:
            indxoaxitemp = retr_indxoaxipnts(gdat, lgal, bgal)
            psfntemp = psfnintp[indxoaxitemp](dist)
        else:
            psfntemp = psfnintp(dist)
    if gdat.kernevaltype == 'bspx':
        pass

    # scale by the PS spectrum
    sbrtpnts = spec[:, None, None] * psfntemp
                
    return sbrtpnts


### find the set of pixels in proximity to a position on the map
def retr_indxpixlevalconc(gdat, strgmodl, dicteval, l, ll, dd):
    
    #print 'll'
    #print ll
    #print 'dd'
    #print dd
    #print 'dicteval[ll][dd][lgal]'
    #summgene(dicteval[ll][dd]['lgal'])
    #print 'dicteval[ll][dd][bgal]'
    #summgene(dicteval[ll][dd]['bgal'])
    #print

    elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    boolelemspec = getattr(gdat, strgmodl + 'boolelemspec')
    
    lgal = dicteval[ll][dd]['lgal']
    bgal = dicteval[ll][dd]['bgal']
    if boolelemspec[l]:
        varbeval = abs(dicteval[ll][dd]['spec'][gdat.indxenerpivt, :])
    if elemtype[l] == 'lens':
        varbeval = dicteval[ll][dd]['defs']
    if elemtype[l] == 'clus':
        varbeval = dicteval[ll][dd]['nobj']
    
    if elemspatevaltype[l] == 'locl' or elemspatevaltype[l] == 'loclhash':
        listindxpixleval = [[] for k in range(lgal.size)]
        for k in range(lgal.size):
            indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
            indxfluxproxtemp = digitize(varbeval[k], gdat.binsprox) - 1
            if gdat.verbtype > 1:
                print 'k'
                print k
                print 'indxpixlpnts'
                print indxpixlpnts
                print 'indxfluxproxtemp'
                print indxfluxproxtemp
            indxpixleval = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
            if isinstance(indxpixleval, int):
                indxpixleval = gdat.indxpixl
            listindxpixleval[k] = indxpixleval
            if gdat.verbtype > 1:
                print 'dicteval[ll][dd][lgal][k]'
                print dicteval[ll][dd]['lgal'][k]
                print 'dicteval[ll][dd][bgal][k]'
                print dicteval[ll][dd]['bgal'][k]
                print 'varbeval'
                print varbeval
                print 'indxpixleval'
                summgene(indxpixleval)
                print
    
    listindxpixlevalconc = unique(concatenate(listindxpixleval))
    
    return listindxpixleval, listindxpixlevalconc


### find the distance between two points on the map
def retr_angldistunit(gdat, lgal, bgal, indxpixleval, retranglcosi=False):
   
    if gdat.pixltype == 'heal':
        xdat, ydat, zaxi = retr_unit(lgal, bgal)
        anglcosi = gdat.xdatgrid[indxpixleval] * xdat + gdat.ydatgrid[indxpixleval] * ydat + gdat.zaxigrid[indxpixleval] * zaxi
        
        if retranglcosi:
            return anglcosi
        else:
            angldist = arccos(anglcosi)
            return angldist
    
    else:
        angldist = sqrt((lgal - gdat.lgalgrid[indxpixleval])**2 + (bgal - gdat.bgalgrid[indxpixleval])**2)
        
        return angldist
    

### find the pixel index of a point on the map
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
    
    # convert to an index of non-zero exposure pixels
    #indxpixl = gdat.indxpixlroficnvt[indxpixl]

    return indxpixl


## obtain count maps
def retr_cntp(gdat, sbrt, indxregieval, indxcubeeval):
   
    cntp = [[] for d in indxregieval]
    for dd, d in enumerate(indxregieval):
        cntp[dd] = sbrt[dd] * gdat.expo[d][indxcubeeval[0][dd]] * gdat.apix
        if gdat.enerdiff:
            cntp[dd] *= gdat.deltener[indxcubeeval[0][dd][0]]
        
    return cntp


def retr_elpsfrac(elpsaxis):
    
    distnorm = sum(((listsamp - gdat.elpscntr[None, :]) / elpsaxis[None, :])**2, axis=1)
    indxsampregu = where(distnorm < 1.)[0]
    thissampfrac = indxsampregu.size / gdat.numbsamp
    vari = (thissampfrac / 0.05 - 1.)**2
    
    return vari


## plotting
### construct path for plots
def retr_plotpath(gdat, gdatmodi, strgstat, strgmodl, strgplot, nameinte=''):
    
    if strgmodl == 'true' or strgstat == '':
        path = gdat.pathinit + nameinte + strgplot + '.pdf'
    elif strgstat == 'post' or strgstat == 'mlik':
        path = gdat.pathplot + gdat.namesampdist + '/finl/' + nameinte + strgstat + strgplot + '.pdf'
    elif strgstat == 'this':
        path = gdat.pathplot + gdat.namesampdist + '/fram/' + nameinte + strgstat + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


### determine the marker size
def retr_mrkrsize(gdat, compampl, namefeatampl):
    
    minm = getattr(gdat, 'minm' + namefeatampl) 
    maxm = getattr(gdat, 'maxm' + namefeatampl)
    mrkrsize = (sqrt(compampl) - sqrt(minm)) / (sqrt(maxm) - sqrt(minm)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
    
    return mrkrsize


## experiment specific

def retr_psfnhubb(gdat):

    # temp
    gdat.psfpexpr = array([0.080, 0.087]) / gdat.anglfact
    gdat.exproaxitype = False


def retr_psfnsdss(gdat):
   
    gdat.psfpexpr = array([0.25 / gdat.anglfact, 1.7e6, 1.9, 0.25 / gdat.anglfact, 2.1e6, 2.])
    gdat.exproaxitype = False


def retr_psfnchan(gdat):

    # temp
    #gdat.psfpexpr = array([0.25, 0.3, 0.4, 0.6, 0.7]) / gdat.anglfact
    if gdat.numbenerfull == 5:
        gdat.psfpexpr = array([0.424 / gdat.anglfact, 2.75, 0.424 / gdat.anglfact, 2.59, 0.440 / gdat.anglfact, 2.47, 0.457 / gdat.anglfact, 2.45, 0.529 / gdat.anglfact, 3.72])
    if gdat.numbenerfull == 2:
        gdat.psfpexpr = array([0.427 / gdat.anglfact, 2.57, 0.449 / gdat.anglfact, 2.49])
    gdat.psfpexpr = gdat.psfpexpr[(2 * gdat.indxenerincl[:, None] + arange(2)[None, :]).flatten()] 
    gdat.exproaxitype = False
    #gdat.psfpexpr = array([0.25 / gdat.anglfact, 
    #                       0.30 / gdat.anglfacti\
    #                       0.40 / gdat.anglfacti\
    #                       0.60 / gdat.anglfacti\
    #                       0.70 / gdat.anglfacti
    #gdat.psfpexpr = array([0.35 / gdat.anglfact, 2e-1, 1.9, 0.5 / gdat.anglfact, 1.e-1, 2.])
    #gdat.psfpexpr = array([0.25 / gdat.anglfact, 2.0e-1, 1.9, \
    #                       0.30 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.40 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.60 / gdat.anglfact, 1.0e-1, 2.0, \
    #                       0.70 / gdat.anglfact, 1.0e-1, 2.0])
    #gdat.exproaxitype = True
   

def retr_psfnsdyn(gdat):

    gdat.psfpexpr = array([0.15 / gdat.anglfact])
    gdat.exproaxitype = False
   

def retr_psfnferm(gdat):
   
    gdat.exproaxitype = False
    
    gdat.recotype = '0007'
    
    if gdat.recotype == '0008':
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
    
    strgpara = ['score', 'gcore', 'stail', 'gtail', 'ntail']
    for m in gdat.indxevtt:
        if gdat.recotype == '0008':
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
            fermform[:, m, k] = interp1d_pick(enerirfn, mean(irfn[strgpara[k]].squeeze(), axis=0))(gdat.meanener)
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 4] = 1. / (1. + fermform[i, m, 4] * fermform[i, m, 2]**2 / fermform[i, m, 0]**2)

    # calculate the scale factor
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # store the fermi PSF parameters
    gdat.psfpexpr = zeros(gdat.numbener * gdat.numbevtt * numbpsfpform)
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            #if k == 0 or k == 2:
            #    gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k] * gdat.fermscalfact[:, m]
            #else:
            #    gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k]
            gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k]
    

def retr_refrchaninit(gdat):
    
    gdat.indxrefr = arange(gdat.numbrefr)
    
    gdat.dictrefr = []
    for q in gdat.indxrefr:
        gdat.dictrefr.append(dict())
    
    gdat.namefeatsignrefr = 'flux'
    
    gdat.legdrefr = ['Xue+2011', 'Wolf+2008']
    
    gdat.listnamefeatamplrefr[0] = 'flux'
    gdat.listnamefeatamplrefr[1] = 'magt'
    gdat.listnamerefr += ['xu11', 'wo08']
    
    setattr(gdat, 'lablotyp', 'O')
    setattr(gdat, 'minmotyp', 0.)
    setattr(gdat, 'maxmotyp', 1.)
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
       
    gdat.refrliststrgfeat[0] += ['lgal', 'bgal', 'flux', 'sind', 'otyp']
    gdat.refrliststrgfeat[1] += ['lgal', 'bgal', 'magt', 'reds', 'otyp']

# col1: Sequential number adopted in this work (Hsu et al. 2014)
# col2-4: ID, R.A.(J2000) and Dec.(J2000) from the CANDELS catalog (Guo et al. 2013).
# col5-7: ID,  R.A.(J2000) and Dec.(J2000) from the MUSYC catalog (Cardamone et al. 2010).
# col8-10: ID, R.A.(J2000) and Dec.(J2000) from the TENIS catalog (Hsieh et al. 2012).
# col11-13: ID, R.A.(J2000) and Dec.(J2000) from the SIMPLE catalog (Damen et al. 2011). 
# col14-17: ID, R.A.(J2000), Dec.(J2000) and positional error from the R13 4Ms-CDFS catalog (Rangel et al. 2013).
# col18-21: ID, R.A.(J2000), Dec.(J2000) and positional error from the X11 4Ms-CDFS catalog (Xue et al. 2011).
# col22-25: ID, R.A.(J2000), Dec.(J2000) and positional error from the L05 250ks-ECDFS catalog (Lehmer et al. 2005).
# col26-29: ID, R.A.(J2000), Dec.(J2000) and positional error from the V06 250ks-ECDFS catalog (Virani et al. 2006).
# col30: "1" indicates that the source is the only possible counterpart to an X-ray source;
#        "n" (2 or more) indicates that the source is one of the n possible counterparts to a X-ray source;
#        "-99" indicates that the source is not associated to any Xray source.
# col31: Posterior value which indicates the reliability of the association to a X-ray source.
#-----------------------------------------------------------------------------------
# [HSN2014]   CANDELS_ID  CANDELS_RA   CANDELS_DEC   MUSYC_ID   MUSYC_RA    MUSYC_DEC    TENIS_ID   TENIS_RA    TENIS_DEC    SIMPLE_ID   SIMPLE_RA   SIMPLE_DEC   R13_ID R13_RA      R13_DEC      R13_PosErr   X11_ID X11_RA      X11_DEC      X11_PosErr   L05_ID L05_RA   L05_DEC   L05_PosErr   V06_ID V06_RA   V06_DEC   V06_PosErr   xflag  post   


def retr_refrchanfinl(gdat):
   
    if gdat.numbsidecart == 300:
        if gdat.anlytype.startswith('extr'):
            gdat.numbpixllgalshft[0] = 1490
            gdat.numbpixlbgalshft[0] = 1430
        if gdat.anlytype.startswith('home'):
            gdat.numbpixllgalshft[0] = 150
            gdat.numbpixlbgalshft[0] = 150
    else:
        raise Exception('Reference elements cannot be aligned with the spatial axes!')
    
    ## WCS object for rotating reference elements into the ROI
    if gdat.numbener == 2:
        gdat.listpathwcss[0] = gdat.pathinpt + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
    else:
        gdat.listpathwcss[0] = gdat.pathinpt + '0.5-0.91028_thresh.img'
    
    # Xue et al. (2011)
    with open(gdat.pathinpt + 'chancatl.txt', 'r') as thisfile:
        pathfile = gdat.pathinpt + 'Xue2011.fits'
        hdun = pf.open(pathfile)
        hdun[0].data
        lgalchan = hdun[1].data['_Glon'] / 180. * pi
        bgalchan = hdun[1].data['_Glat'] / 180. * pi
        fluxchansoft = hdun[1].data['SFlux']
        fluxchanhard = hdun[1].data['HFlux']
        objttypechan = hdun[1].data['Otype']
    
    # position
    gdat.refrlgal[0][0] = lgalchan
    gdat.refrbgal[0][0] = bgalchan

    # spectra
    gdat.refrspec = [[zeros((3, gdat.numbener, lgalchan.size))]]
    if gdat.numbener == 2:
        gdat.refrspec[0][0][0, 0, :] = fluxchansoft * 0.624e9
        gdat.refrspec[0][0][0, 1, :] = fluxchanhard * 0.624e9 / 16.
    else:
        gdat.refrspec[0][0][0, :, :] = fluxchansoft[None, :] * 0.624e9
    gdat.refrspec[0][0][1, :, :] = gdat.refrspec[0][0][0, :, :]
    gdat.refrspec[0][0][2, :, :] = gdat.refrspec[0][0][0, :, :]
   
    # fluxes
    gdat.refrflux[0][0] = gdat.refrspec[0][0][:, gdat.indxenerpivt, :]

    # spectral indices
    if gdat.numbener > 1:
        gdat.refrsind[0][0] = -log(gdat.refrspec[0][0][0, 1, :] / gdat.refrspec[0][0][0, 0, :]) / log(sqrt(7. / 2.) / sqrt(0.5 * 2.))
        gdat.refrsind[0][0][where(logical_not(isfinite(gdat.refrsind[0][0])))[0]] = 2.

    ## object type
    indx = where(objttypechan == 'AGN')[0]
    objttypechan[indx] = 0
    indx = where(objttypechan == 'Galaxy')[0]
    objttypechan[indx] = 1
    indx = where(objttypechan == 'Star')[0]
    objttypechan[indx] = 2
    objttypechan = objttypechan.astype(int)
    gdat.refrotyp[0][0] = objttypechan.astype(float)

    # Wolf et al. (2011)
    path = gdat.pathdata + 'inpt/Wolf2008.fits'
    data = pf.getdata(path)
    gdat.refrlgal[1][0] = deg2rad(data['_Glon'])
    gdat.refrlgal[1][0] = ((gdat.refrlgal[1][0] - pi) % (2. * pi)) - pi
    gdat.refrbgal[1][0] = deg2rad(data['_Glat'])
    gdat.refrmagt[1][0] = data['Rmag']
    gdat.refrreds[1][0] = data['MCz']
    #listname = []
    #for k in range(data['MCclass'].size):
    #    if not data['MCclass'][k] in listname:
    #        listname.append(data['MCclass'][k])
    listname = ['Galaxy', 'Galaxy  (Uncl!)', 'QSO     (Gal?)', 'Galaxy  (Star?)', 'Star', 'Strange Object', 'QSO', 'WDwarf']
    gdat.refrotyp[1][0] = empty_like(gdat.refrreds[1][0]) 
    for k, name in enumerate(listname):
        indx = where(gdat.refrotyp[1][0] == name)
        gdat.refrotyp[1][0][indx] = k
     

def retr_refrferminit(gdat):
    
    gdat.listnamerefr += ['ac15', 'ma05']
    gdat.indxrefr = arange(gdat.numbrefr)
    
    gdat.legdrefr = ['Acero+2015', 'Manchester+2005']

    gdat.listnamefeatamplrefr[0] = 'flux'
    gdat.listnamefeatamplrefr[1] = 'flux0400'
    gdat.namefeatsignrefr = 'flux'
    
    setattr(gdat, 'lablcurvac15', '%s_{3FGL}' % gdat.lablcurv)
    setattr(gdat, 'lablexpcac15', 'E_{c,3FGL}')
    
    for name in gdat.listnamerefr:
        setattr(gdat, 'minmcurv' + name, -1.)
        setattr(gdat, 'maxmcurv' + name, 1.)
        setattr(gdat, 'minmexpc' + name, 0.1)
        setattr(gdat, 'maxmexpc' + name, 10.)
   
    gdat.refrliststrgfeat[0] += ['lgal', 'bgal', 'flux', 'sind', 'curv', 'expc']
    gdat.refrliststrgfeat[1] += ['lgal', 'bgal', 'flux0400', 'per0', 'per1']


def retr_refrfermfinl(gdat):
    
    # Acero+2015
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'
    fgl3 = pf.getdata(path)
    
    gdat.refrlgal[0][0] = deg2rad(fgl3['glon'])
    gdat.refrlgal[0][0] = ((gdat.refrlgal[0][0] - pi) % (2. * pi)) - pi
    gdat.refrbgal[0][0] = deg2rad(fgl3['glat'])
    
    gdat.refrnumbelemfull = gdat.refrlgal[0][0].size

    gdat.refrspec = [[empty((3, gdat.numbener, gdat.refrlgal[0][0].size))]]
    gdat.refrspec[0][0][0, :, :] = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], \
                                                                                            fgl3['Flux10000_100000']))[gdat.indxenerincl, :] / gdat.deltener[:, None]
    
    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.refrspec[0][0][1, :, :] = gdat.refrspec[0][0][0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.refrspec[0][0][2, :, :] = gdat.refrspec[0][0][0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.refrspec[0][0][where(isfinite(gdat.refrspec[0][0]) == False)] = 0.
    
    gdat.refrflux[0][0] = gdat.refrspec[0][0][:, gdat.indxenerpivt, :]
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3strg = fgl3['Source_Name']
    fgl3strgclss = fgl3['CLASS1']
    fgl3strgassc = fgl3['ASSOC1']
    
    fgl3timevari = fgl3['Variability_Index']

    fgl3spectype = fgl3['SpectrumType']
    gdat.refrsind[0][0] = fgl3['Spectral_Index']
    gdat.refrcurv[0][0] = fgl3['beta']
    gdat.refrexpc[0][0] = fgl3['Cutoff'] * 1e-3
    
    gdat.refrcurv[0][0][where(logical_not(isfinite(gdat.refrcurv[0][0])))] = -10.
    gdat.refrexpc[0][0][where(logical_not(isfinite(gdat.refrexpc[0][0])))] = 0.
    
    gdat.refrsind[0][0] = tile(gdat.refrsind[0][0], (3, 1)) 
    gdat.refrcurv[0][0] = tile(gdat.refrcurv[0][0], (3, 1)) 
    gdat.refrexpc[0][0] = tile(gdat.refrexpc[0][0], (3, 1)) 

    #indxtimevari = where(fgl3timevari < 100.)[0]
    #indxtimevari = where(gdat.refrspec[0][0, gdat.indxenerpivt, :] > gdat.minmflux)[0]

    # Manchester+2005
    path = gdat.pathdata + 'inpt/Manchester2005.fits'
    data = pf.getdata(path)
   
    gdat.refrlgal[1][0] = deg2rad(data['glon'])
    gdat.refrlgal[1][0] = ((gdat.refrlgal[1][0] - pi) % (2. * pi)) - pi
    gdat.refrbgal[1][0] = deg2rad(data['glat'])
    
    gdat.refrper0[1][0] = data['P0']
    gdat.refrper1[1][0] = data['P1']
    gdat.refrflux0400[1][0] = data['S400']
    #gdat.refrdism[1][0] = data['DM']
    #gdat.refrdlos[1][0] = data['Dist']

    setattr(gdat, 'lablper0', 'P_0')
    setattr(gdat, 'minmper0', 1e-4)
    setattr(gdat, 'maxmper0', 1e1)
    setattr(gdat, 'factper0plot', 1.)
    setattr(gdat, 'scalper0plot', 'logt')
    
    setattr(gdat, 'lablper1', 'P_1')
    setattr(gdat, 'minmper1', 1e-20)
    setattr(gdat, 'maxmper1', 1e-10)
    setattr(gdat, 'factper1plot', 1.)
    setattr(gdat, 'scalper1plot', 'logt')
    
    setattr(gdat, 'lablflux0400', 'S_{400}')
    setattr(gdat, 'minmflux0400', 1e-1)
    setattr(gdat, 'maxmflux0400', 1e4)
    setattr(gdat, 'factflux0400plot', 1.)
    setattr(gdat, 'scalflux0400plot', 'logt')
    
    # error budget
    for name in ['lgal', 'bgal', 'per0', 'per1', 'flux0400']:
        refrtile = [[[] for d in gdat.indxregi] for q in gdat.indxrefr]
        refrfeat = getattr(gdat, 'refr' + name)
        for q in gdat.indxrefr:
            for d in gdat.indxregi:
                if len(refrfeat[q][d]) > 0:
                    refrtile[q][d] = tile(refrfeat[q][d], (3, 1))
        setattr(gdat, 'refr' + name, refrtile)


### find the list of pairs for splits and merges
def retr_listpair(gdat, lgal, bgal):
    
    if gdat.verbtype > 1:
        print 'Finding element pairs inside the linking length...'
    
    listpair = []
    for k in range(lgal.size):
        # temp -- linking uses the Cartesian approximation, which is accurate enough for splits and merges inside a small circle
        #indxelem = k + 1 + where(sqrt((bgal[k+1:] - bgal[k])**2 + (lgal[k+1:] - lgal[k])**2) < gdat.radispmr)[0]
        indxelem = k + 1 + where((fabs(bgal[k+1:] - bgal[k]) < gdat.radispmr) & (fabs(lgal[k+1:] - lgal[k]) < gdat.radispmr))[0]
        for n in range(indxelem.size):
            listpair.append([k, indxelem[n]])
    
    if gdat.verbtype > 1:
        numbpair = len(listpair)
        print '%d pairs found.' % numbpair

    if gdat.diagmode:
        boolgood = True
        for n in range(len(listpair)):
            dist = sqrt((lgal[listpair[n][0]] - lgal[listpair[n][1]])**2 + (bgal[listpair[n][0]] - bgal[listpair[n][1]])**2)
            if dist >= gdat.radispmr:
                boolgood = False
        if not boolgood:
            Exception('Inappropriate list of pairs')

    return listpair


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


def show_samp(gdat, gdatmodi, thisonly=False):
    
    if not thisonly:
        print '%5s %22s %14s %14s %14s %14s %14s %14s %10s' % ('index', 'name', 'thissamp', 'nextsamp', 'thissampvarb', 'nextsampvarb', 'diffsampvarb', 'prop', 'indxstdp')
    else:
        print '%5s %22s %14s %14s' % ('index', 'name', 'thissamp', 'thissampvarb')
    
    for k in gdat.fittindxpara:
        if k == gdat.fittnumbfixp:
            print
        if k < gdat.fittnumbfixp:
            name = gdat.fittnamefixp[k]
        else:
            name = ''

        if gdat.fittnumbtrap > 0:
            gdatmodi.indxparaprop = concatenate((gdat.indxfixpprop, concatenate(concatenate(gdatmodi.thisindxsampcomp['comp']))))
        else:
            gdatmodi.indxparaprop = gdat.indxfixpprop
        try:
            strgboolmodi = '%s' % (k in gdatmodi.indxsampmodi)
        except:
            strgboolmodi = ''
        if not thisonly:
            print '%5d %22s %14.4g %14.4g %14.4g %14.4g %14.4g %14s %10d' % (k, name, getattr(gdatmodi, 'thissamp')[k], getattr(gdatmodi, 'nextsamp')[k], gdatmodi.thissampvarb[k], \
                                               gdatmodi.nextsampvarb[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k], strgboolmodi, gdat.indxstdppara[k])
        else:
            print '%5d %22s %14.4g %14.4g' % (k, name, getattr(gdatmodi, 'thissamp')[k], gdatmodi.thissampvarb[k])
    

def show_samplong(gdat, gdatmodi):
    
    for strgmodl in gdat.liststrgmodl:
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
            for d in gdat.fittindxregipopl[l]:
                comp = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[strgcomp][l][d]]
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
                    gdatmodi.nextsamp[gdatmodi.thisindxsampcomp[strgcomp][l][d]] = unit


def prop_stat(gdat, gdatmodi, strgmodl, thisindxelem=None, thisindxpopl=None, thisindxregi=None, brth=False, deth=False):
 
    if gdat.verbtype > 1:
        print 'prop_stat()'
    
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    indxcomp = getattr(gdat, strgmodl + 'indxcomp')
    namepara = getattr(gdat, strgmodl + 'namepara')
    minmnumbelem = getattr(gdat, strgmodl + 'minmnumbelem')
    maxmnumbelem = getattr(gdat, strgmodl + 'maxmnumbelem')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    if gdat.numbpixl > 1:
        indxfixppsfp = getattr(gdat, strgmodl + 'indxfixppsfp')
    indxfixpbacp = getattr(gdat, strgmodl + 'indxfixpbacp')
    indxfixplenp = getattr(gdat, strgmodl + 'indxfixplenp')
    indxfixpsour = getattr(gdat, strgmodl + 'indxfixpsour')
    indxfixphost = getattr(gdat, strgmodl + 'indxfixphost')
    if numbtrap > 0:
        indxfixpdist = getattr(gdat, strgmodl + 'indxfixpdist')
        indxfixpmeanelem = getattr(gdat, strgmodl + 'indxfixpmeanelem')
        indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
        boolelempsfn = getattr(gdat, strgmodl + 'boolelempsfn')
    indxpopl = getattr(gdat, strgmodl + 'indxpopl')
    indxregipopl = getattr(gdat, strgmodl + 'indxregipopl')
    namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
    indxsampcompinit = getattr(gdat, strgmodl + 'indxsampcompinit')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    numbtrappoplcumr = getattr(gdat, strgmodl + 'numbtrappoplcumr')
    numbtrapregipoplcuml = getattr(gdat, strgmodl + 'numbtrapregipoplcuml')
    numbtrapregipoplcumr = getattr(gdat, strgmodl + 'numbtrapregipoplcumr')
    elemregitype = getattr(gdat, strgmodl + 'elemregitype')
    
    spectype = getattr(gdat, strgmodl + 'spectype')
        
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, 'this', strgmodl)

    strgpfixthis = retr_strgpfix('this', strgmodl)
    strgpfixnext = retr_strgpfix('next', strgmodl)
    
    thissamp = getattr(gdatobjt, strgpfixthis + 'samp')
    thissampvarb = getattr(gdatobjt, strgpfixthis + 'sampvarb')
    if numbtrap > 0:
        boolelemspec = getattr(gdat, strgmodl + 'boolelemspec')
        thisindxelemfull = getattr(gdatobjt, strgpfixthis + 'indxelemfull')
        thisindxsampcomp = retr_indxsampcomp(gdat, thisindxelemfull, strgmodl)
        setattr(gdatobjt, strgpfixthis + 'indxsampcomp', thisindxsampcomp)
    
    nextsamp = copy(thissamp)
    nextsampvarb = copy(thissampvarb)
    if numbtrap > 0:
        nextindxelemfull = deepcopy(thisindxelemfull)
     
    setattr(gdatobjt, strgpfixnext + 'samp', nextsamp)
    setattr(gdatobjt, strgpfixnext + 'sampvarb', nextsampvarb)
    if numbtrap > 0:
        setattr(gdatobjt, strgpfixnext + 'indxelemfull', nextindxelemfull)
    
    gdatmodi.thisaccppsfn = True
    gdatmodi.thisaccpprio = True
    
    gdatmodi.evalllikpert = False 

    # initialize the Boolean flag indicating the type of transdimensional proposal
    gdatmodi.propcomp = False
    gdatmodi.propfixp = False
    gdatmodi.propbacp = False
    gdatmodi.proppsfp = False
    gdatmodi.proplenp = False
    gdatmodi.propwith = False
    gdatmodi.propbrth = False
    gdatmodi.propdeth = False
    gdatmodi.propsplt = False
    gdatmodi.propmerg = False
    gdatmodi.propmeanelem = False
    gdatmodi.propdist = False

    gdatmodi.nextljcbfact = 0.
    gdatmodi.nextlpautotl = 0. 
    gdatmodi.nextlcomfact = 0.
   
    gdatmodi.prophost = False
    gdatmodi.propsour = False
    gdatmodi.proppsfnconv = False
    
    # index of the population in which a transdimensional proposal will be made
    if thisindxpopl == None:
        gdatmodi.indxpoplmodi = array([choice(indxpopl)])
    else:
        gdatmodi.indxpoplmodi = array([thisindxpopl])
    
    if thisindxregi == None:
        if elemregitype[gdatmodi.indxpoplmodi[0]]:
            gdatmodi.indxregimodi = array([choice(gdat.indxregi)])
        else:
            gdatmodi.indxregimodi = array([0])
    else:
        gdatmodi.indxregimodi = array([thisindxregi])
    
    # forced death or birth does not check for the prior on the dimensionality on purpose!
    if deth or brth or rand() < gdat.probtran and ( \
                                                 thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] != \
                                                 minmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]] or \
                                                 thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] != \
                                                 maxmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]):
        
        if brth or deth or rand() < gdat.probbrde / gdat.probtran:
            
            ## births and deaths
            if thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] == \
                                                                                            maxmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]] or deth:
                gdatmodi.propdeth = True
            elif thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] == \
                                                                                            minmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]] or brth:
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
            if thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] == maxmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]:
                gdatmodi.propmerg = True
            elif thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] == minmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]:
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
        gdatmodi.indxenermodi = gdat.indxener
        gdatmodi.indxevttmodi = gdat.indxevtt

    else:
        if gdat.propfixp and gdat.propcomp and gdat.indxfixpprop.size > 0:
            gdatmodi.thisindxsampfull = concatenate((gdat.indxfixpprop, concatenate(concatenate(thisindxsampcomp['comp']))))
        elif gdat.propcomp and not gdat.propfixp:
            gdatmodi.thisindxsampfull = concatenate(concatenate(thisindxsampcomp['comp']))
        else:
            gdatmodi.thisindxsampfull = gdat.indxfixpprop
        gdatmodi.indxsampmodi = choice(gdatmodi.thisindxsampfull)
        
        if gdatmodi.indxsampmodi in indxfixp:
            gdatmodi.propfixp = True
        else:
            gdatmodi.propcomp = True
            
        if gdat.fittpsfnevaltype != 'none' and gdatmodi.indxsampmodi in indxfixppsfp:
            gdatmodi.proppsfp = True
            gdatmodi.indxregimodi = gdat.indxregi
            gdatmodi.indxenermodi = array([int(gdat.fittnamepara[gdatmodi.indxsampmodi][-5])])
            gdatmodi.indxevttmodi = array([int(gdat.fittnamepara[gdatmodi.indxsampmodi][-1])])
            
        elif gdatmodi.indxsampmodi in indxfixpbacp:
            gdatmodi.propbacp = True
            gdatmodi.indxbackmodi = array([int(gdat.fittnamepara[gdatmodi.indxsampmodi][8:12])])
            gdatmodi.indxregimodi = array([int(gdat.fittnamepara[gdatmodi.indxsampmodi][15])])
            if gdat.fittspecback[gdatmodi.indxbackmodi[0]]:
                gdatmodi.indxenermodi = gdat.indxener
            else:
                gdatmodi.indxenermodi = array([int(gdat.fittnamepara[gdatmodi.indxsampmodi][19])])
            gdatmodi.indxevttmodi = gdat.indxevtt
            
        elif gdatmodi.indxsampmodi in indxfixplenp:
            gdatmodi.indxregimodi = array([int(gdat.fittnamepara[gdatmodi.indxsampmodi][-1])])
            gdatmodi.indxenermodi = gdat.indxener
            gdatmodi.indxevttmodi = gdat.indxevtt
            
        if False and (gdatmodi.indxsampmodi in indxfixpbacp or gdat.fittpsfnevaltype != 'none' and gdatmodi.indxsampmodi in indxfixppsfp):
            print 'gdat.fittnamepara[gdatmodi.indxsampmodi]'
            print gdat.fittnamepara[gdatmodi.indxsampmodi]
            print 'gdatmodi.indxenermodi'
            print gdatmodi.indxenermodi
            print 'gdatmodi.indxevttmodi'
            print gdatmodi.indxevttmodi
            print 'gdatmodi.propbacp'
            print gdatmodi.propbacp
            print 'gdatmodi.proppsfp'
            print gdatmodi.proppsfp
                
        # temp
        if gdat.fittnumbtrap > 0:
            if gdatmodi.indxsampmodi in indxfixpmeanelem:
                gdatmodi.propmeanelem = True
            if gdatmodi.indxsampmodi in indxfixpdist:
                gdatmodi.propdist = True
        
        gdatmodi.thisindxproptype = gdat.indxstdppara[gdatmodi.indxsampmodi]
        gdatmodi.propwith = True
        gdatmodi.proppsfnconv = (gdat.fittpsfnevaltype == 'conv' or gdat.fittpsfnevaltype == 'full') and gdatmodi.proppsfp
        if gdat.fittlensmodltype != 'none' or gdat.fitthostemistype != 'none':
            gdatmodi.prophost = gdatmodi.indxsampmodi in indxfixphost
        if gdat.fittlensmodltype != 'none':
            gdatmodi.propsour = gdatmodi.indxsampmodi in indxfixpsour
        if gdat.fittlensmodltype != 'none':
            gdatmodi.indxenermodi = gdat.indxener
            gdatmodi.indxevttmodi = gdat.indxevtt
    
    # derived proposal flags
    gdatmodi.proptran = gdatmodi.propbrth or gdatmodi.propdeth or gdatmodi.propsplt or gdatmodi.propmerg
    gdatmodi.propelem = gdatmodi.propcomp or gdatmodi.proptran
    gdatmodi.prophypr = gdatmodi.propmeanelem or gdatmodi.propdist
    gdatmodi.proplpri = gdatmodi.prophypr
    gdatmodi.propllik = not gdatmodi.prophypr
    gdatmodi.evalllikpert = numbtrap > 0 and (not gdatmodi.propfixp or gdatmodi.proppsfp and True in boolelempsfn)

    if gdat.verbtype > 1:
        if gdat.fittnumbtrap > 0:
            print 'thisindxelemfull'
            print thisindxelemfull
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
        if not gdatmodi.propfixp:
            print 'gdatmodi.indxpoplmodi'
            print gdatmodi.indxpoplmodi
        print 'gdatmodi.evalllikpert'
        print gdatmodi.evalllikpert
        print 'gdatmodi.propfixp'
        print gdatmodi.propfixp
        print 'gdatmodi.propbacp'
        print gdatmodi.propbacp
        print 'gdatmodi.proppsfp'
        print gdatmodi.proppsfp
        print 'gdatmodi.prophost'
        print gdatmodi.prophost
        print 'gdatmodi.propcomp'
        print gdatmodi.propcomp
        if gdatmodi.propwith:
            print 'indxsampmodi'
            print gdatmodi.indxsampmodi
    
    if gdat.diagmode:
        if gdat.probbrde == 0. and (gdatmodi.propbrth or gdatmodi.propdeth and not deth):
            raise Exception('')

        if gdat.probspmr == 0. and (gdatmodi.propsplt or gdatmodi.propmerg):
            raise Exception('')

    stdvstdp = gdat.stdvstdp * gdatmodi.thisstdpscalfact * gdatmodi.thistmprfactstdv

    gdatmodi.numbelemeval = [[[]]]
    if gdatmodi.propwith:
        
        if not gdatmodi.propfixp:
            
            # given the index of the element parameter to be perturbed, find the population, region, element and component indices
            gdatmodi.indxtrapmodi = gdatmodi.indxsampmodi - indxsampcompinit
            gdatmodi.indxpoplmodi = array([amin(where(gdatmodi.indxtrapmodi // numbtrappoplcumr == 0))])
            gdatmodi.indxregimodi = array([amin(where(gdatmodi.indxtrapmodi // numbtrapregipoplcumr[gdatmodi.indxpoplmodi[0]] == 0))])
            gdatmodi.numbparapoplinit = gdatmodi.indxtrapmodi - numbtrapregipoplcuml[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]
            gdatmodi.indxelemmodi = [gdatmodi.numbparapoplinit // numbcomp[gdatmodi.indxpoplmodi[0]]]
            gdatmodi.indxcompmodi = gdatmodi.numbparapoplinit % numbcomp[gdatmodi.indxpoplmodi[0]]
            
            # indices of the energy and PSF bins over which likelihood will be evaluated
            gdatmodi.indxenermodi = gdat.indxener
            gdatmodi.indxevttmodi = gdat.indxevtt
            
            # number of evaluation elements
            gdatmodi.numbelemeval[0][0] = 2
            
            gdatmodi.indxelemfullmodi = [gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]].index(gdatmodi.indxelemmodi[0])]
            gdatmodi.thisstrgcomp = liststrgcomp[gdatmodi.indxpoplmodi[0]][gdatmodi.indxcompmodi]
            gdatmodi.thisliststrgcomp = [[] for l in indxpopl]
            gdatmodi.thisliststrgcomp[gdatmodi.indxpoplmodi[0]].append(gdatmodi.thisstrgcomp)
            gdatmodi.thislistscalcomp = [[] for l in indxpopl]
            gdatmodi.thislistscalcomp[gdatmodi.indxpoplmodi[0]].append(listscalcomp[gdatmodi.indxpoplmodi[0]][gdatmodi.indxcompmodi])
            gdatmodi.thisamplpert = thissampvarb[thisindxsampcomp[gdat.fittnamefeatampl[gdatmodi.indxpoplmodi[0]]][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]], None]
        else:
            gdatmodi.numbelemeval[0][0] = 0
            
        ## index of the proposal
        gdatmodi.indxstdpmodi = gdat.indxstdppara[gdatmodi.indxsampmodi]
        ## proposal scale
        if gdatmodi.propfixp:
            gdatmodi.thisstdvpara = stdvstdp[gdatmodi.indxstdpmodi]
        else:
            stdvstdpcomp = stdvstdp[gdatmodi.indxstdpmodi]
            if False:
                # parameter-dependent proposal scale
                thiscompampl = thissampvarb[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][d][indxelemfull]]
                nextcompampl = gdatmodi.nextsampvarb[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][d][indxelemfull]]
                minmcompampl = getattr(gdat, 'minm' + gdat.fittnamefeatampl[l])
                thiscompunit = gdatmodi.thissamp[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][d][indxelemfull]]
                nextcompunit = gdatmodi.nextsamp[thisindxsampcomp[gdat.fittnamefeatampl[l]][l][d][indxelemfull]]
                if strgcomp == gdat.fittnamefeatampl[l]:
                    # temp -- this only works if compampl is powr distributed
                    gdatmodi.thisstdvpara = stdvstdpcomp / (thiscompampl / minmcompampl)**2.
                    gdatmodi.nextstdv = stdvstdpcomp / (nextcompampl / minmcompampl)**2.
                    gdatmodi.thislfctasym += sum(0.5 * (nextcompunit - thiscompunit)**2 * (1. / gdatmodi.thisstdv**2 - 1. / gdatmodi.nextstdv**2))
                else:
                    gdatmodi.thisstdvpara = stdvstdpcomp / (minimum(thiscompampl, nextcompampl) / minmcompampl)**0.5
            else:
                gdatmodi.thisstdvpara = stdvstdpcomp
        ## propose a step
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdatmodi.thisstdvpara)
        
        if gdatmodi.propfixp:
            gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_fixp(gdat, strgmodl, gdatmodi.nextsamp[gdatmodi.indxsampmodi], gdatmodi.indxsampmodi)
            
        # temp
        if numbtrap > 0 and gdatmodi.propdist:
            # temp
            # this should check whether rescaling is necessary
            gdatmodi.thisstrgcomp = namepara[gdatmodi.indxsampmodi][-16:-12]

            ### rescale the element components due to the hyperparameter step
            rscl_elem(gdat, gdatmodi, gdatmodi.indxsampmodi)
        
        if gdatmodi.propcomp:
            retr_sampvarbtrap(gdat, strgmodl, thisindxsampcomp, gdatmodi.indxregimodi, gdatmodi.indxpoplmodi, gdatmodi.thisliststrgcomp, \
                                                          gdatmodi.thislistscalcomp, gdatmodi.nextsamp, gdatmodi.nextsampvarb, indxelemfull=gdatmodi.indxelemfullmodi[0])
            
            gdatmodi.thiscomp = thissampvarb[thisindxsampcomp[gdatmodi.thisstrgcomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
            gdatmodi.nextcomp = nextsampvarb[thisindxsampcomp[gdatmodi.thisstrgcomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
        

    else:
        # number of unit sample vector elements to be modified
        numbcompmodi = numbcomp[gdatmodi.indxpoplmodi[0]]
   
    if gdat.verbtype > 1:
        if gdatmodi.propllik:
            print 'gdatmodi.indxregimodi'
            print gdatmodi.indxregimodi
            print 'gdatmodi.indxenermodi'
            print gdatmodi.indxenermodi
            print 'gdatmodi.indxevttmodi'
            print gdatmodi.indxevttmodi
    
    if gdatmodi.proptran:
        gdatmodi.auxipara = rand(numbcomp[gdatmodi.indxpoplmodi[0]])
        gdatmodi.indxsamptran = []
    
    if gdatmodi.propbrth or gdatmodi.propsplt:
       
        # find an empty slot in the element list
        for k in range(maxmnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]):
            if not k in gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]:
                break
       
        # sample indices to add the new element
        gdatmodi.indxsamptran.append(retr_indxsamppnts(gdat, gdatmodi.indxpoplmodi[0], gdatmodi.indxregimodi[0], array([k])))
        gdatmodi.indxelemmodi = [k]
        gdatmodi.nextindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]].append(k)
        gdatmodi.indxelemfullmodi = [thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]].astype(int)]
    if gdatmodi.propbrth:
        
        gdatmodi.numbelemeval[0][0] = 1
        
        # sample auxiliary variables
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.auxipara
    
    # death
    if gdatmodi.propdeth:
        
        gdatmodi.numbelemeval[0][0] = 1
        
        # occupied element index to be killed
        if thisindxelem == None:
            dethindxindxelem = choice(arange(thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]], dtype=int))
        else:
            dethindxindxelem = thisindxelem

        # element index to be killed
        gdatmodi.indxelemmodi = []
        gdatmodi.indxelemfullmodi = []
        if gdat.verbtype > 1:
            print 'dethindxindxelem'
            print dethindxindxelem

        gdatmodi.indxelemmodi.append(thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][dethindxindxelem])
        gdatmodi.indxelemfullmodi.append(dethindxindxelem)
        # parameter indices to be killed
        gdatmodi.indxsamptran.append(indxsampcompinit + numbtrapregipoplcuml[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]] + \
                                                    gdatmodi.indxelemmodi[0] * numbcomp[gdatmodi.indxpoplmodi[0]] + indxcomp[gdatmodi.indxpoplmodi[0]])
    
    if gdatmodi.propsplt or gdatmodi.propmerg:
        gdatmodi.comppare = empty(numbcomp[gdatmodi.indxpoplmodi[0]])
        gdatmodi.compfrst = empty(numbcomp[gdatmodi.indxpoplmodi[0]])
        gdatmodi.compseco = empty(numbcomp[gdatmodi.indxpoplmodi[0]])
    
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbelemeval[0][0] = 3
        
        gdatmodi.indxelemfullsplt = choice(arange(gdatmodi.thissamp[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]], dtype=int))
        gdatmodi.indxelemsplt = thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullsplt]
        gdatmodi.indxelemfullmodi.insert(0, gdatmodi.indxelemfullsplt)
        gdatmodi.indxelemmodi.insert(0, gdatmodi.indxelemsplt)

        # boundaries of the sample vector indices to be modified
        ## first
        gdatmodi.indxsampfrst = indxsampcompinit + numbtrap * gdatmodi.indxpoplmodi[0] + gdatmodi.indxelemsplt * numbcomp[gdatmodi.indxpoplmodi[0]]
        gdatmodi.indxsampfrstfinl = gdatmodi.indxsampfrst + numbcomp[gdatmodi.indxpoplmodi[0]]
        gdatmodi.indxsamptran.insert(0, arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl))
        
        ## second
        gdatmodi.indxsampseco = indxsampcompinit + numbtrap * gdatmodi.indxpoplmodi[0] + gdatmodi.indxelemmodi[1] * numbcomp[gdatmodi.indxpoplmodi[0]]
        gdatmodi.indxsampsecofinl = gdatmodi.indxsampseco + numbcomp[gdatmodi.indxpoplmodi[0]]
        
        # concatenated sample vector indices to be modified
        gdatmodi.indxsampmodi = concatenate((arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl, dtype=int), \
                                                                                    arange(gdatmodi.indxsampseco, gdatmodi.indxsampsecofinl, dtype=int)))
        
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            gdatmodi.comppare[k] = \
                    copy(thissampvarb[thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]])

        # determine the new element parameters
        # temp -- only valid for power-law energy spectrum
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            if k == 0 or k == 1:
                gdatmodi.auxipara[k] = randn() * gdat.radispmr
            elif k == 2:
                gdatmodi.auxipara[k] = rand()
            else:
                if listscalcomp[gdatmodi.indxpoplmodi[0]] == 'self':
                    minm = getattr(gdat, strgmodl + 'minm' + strgcomp)
                    fact = getattr(gdat, strgmodl + 'fact' + strgcomp)
                    gdatmodi.auxipara[k] = icdf_self(rand(), minm, fact)
                elif listscalcomp[gdatmodi.indxpoplmodi[0]] == 'gaus':
                    distmean = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'distmean')[gdatmodi.indxpoplmodi[0]]]
                    diststdv = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgcomp + 'diststdv')[gdatmodi.indxpoplmodi[0]]]
                    gdatmodi.auxipara[k] = icdf_gaus(rand(), distmean, diststdv)
        
        gdatmodi.compfrst[0] = gdatmodi.comppare[0] + (1. - gdatmodi.auxipara[2]) * gdatmodi.auxipara[0]
        gdatmodi.compfrst[1] = gdatmodi.comppare[1] + (1. - gdatmodi.auxipara[2]) * gdatmodi.auxipara[1]
        gdatmodi.compfrst[2] = gdatmodi.auxipara[2] * gdatmodi.comppare[2]
        gdatmodi.compseco[0] = gdatmodi.comppare[0] - gdatmodi.auxipara[2] * gdatmodi.auxipara[0]
        gdatmodi.compseco[1] = gdatmodi.comppare[1] - gdatmodi.auxipara[2] * gdatmodi.auxipara[1]
        gdatmodi.compseco[2] = (1. - gdatmodi.auxipara[2]) * gdatmodi.comppare[2]
        #for k in range(numbcomp[gdatmodi.indxpoplmodi[0]]):
        #    if k > 2:
        #        gdatmodi.compfrst[k] = gdatmodi.comppare[k]
        #        gdatmodi.compseco[k] = gdatmodi.comppare[k]
        
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0][:3]] = cdfn_trap(gdat, gdatmodi, gdatmodi.compfrst)[:3]
        gdatmodi.nextsamp[gdatmodi.indxsamptran[1][:3]] = cdfn_trap(gdat, gdatmodi, gdatmodi.compseco)[:3]
        gdatmodi.nextsampvarb[gdatmodi.indxsamptran[0][:3]] = gdatmodi.compfrst[:3]
        gdatmodi.nextsampvarb[gdatmodi.indxsamptran[1][:3]] = gdatmodi.compseco[:3]
        
        if fabs(gdatmodi.compfrst[0]) > gdat.maxmgang or fabs(gdatmodi.compseco[0]) > gdat.maxmgang or \
                     fabs(gdatmodi.compfrst[1]) > gdat.maxmgang or fabs(gdatmodi.compseco[1]) > gdat.maxmgang or \
                     gdatmodi.compfrst[2] < getattr(gdat, strgmodl + 'minm' + gdat.fittnamefeatampl[gdatmodi.indxpoplmodi[0]]) or \
                     gdatmodi.compseco[2] < getattr(gdat, strgmodl + 'minm' + gdat.fittnamefeatampl[gdatmodi.indxpoplmodi[0]]):
            gdatmodi.thisaccpprio = False

        if False or gdat.verbtype > 1:
            print 'gdatmodi.indxsamptran'
            print gdatmodi.indxsamptran
            print 'gdatmodi.auxipara'
            print gdatmodi.auxipara
            for k in range(3):
                print 'k'
                print k
                print 'gdatmodi.compfrst[k]'
                print gdatmodi.compfrst[k]
                print 'gdatmodi.compseco[k]'
                print gdatmodi.compseco[k]
                print 'gdatmodi.comppare[k]'
                print gdatmodi.comppare[k]
                print
            print 'gdatmodi.compfrst[2] + gdatmodi.compseco[2]'
            print gdatmodi.compfrst[2] + gdatmodi.compseco[2]
            if abs(gdatmodi.compfrst[2] + gdatmodi.compseco[2] - 3e-8) / (gdatmodi.compfrst[2] + gdatmodi.compseco[2]) > 0.001:
                raise Exception('')
            if abs(gdatmodi.comppare[2] - 3e-8) / gdatmodi.comppare[2] > 0.001:
                raise Exception('')
            if 100. * (gdatmodi.compfrst[2] + gdatmodi.compseco[2] - gdatmodi.comppare[2]) / gdatmodi.comppare[2] > 1.:
                print 'gdatmodi.auxipara[2]'
                print gdatmodi.auxipara[2]
                print 'gdatmodi.comppare[2]'
                print gdatmodi.comppare[2]
                raise Exception('')
            print

        # calculate the list of pairs
        if False and gdatmodi.thisaccpprio:

            # calculate the list of pairs
            ## proposed
            lgal = concatenate((array([gdatmodi.compfrst[0], gdatmodi.compseco[0]]), \
                           setdiff1d(thissampvarb[thisindxsampcomp['lgal'][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]], gdatmodi.comppare[0])))
            bgal = concatenate((array([gdatmodi.compfrst[1], gdatmodi.compseco[1]]), \
                           setdiff1d(thissampvarb[thisindxsampcomp['bgal'][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]], gdatmodi.comppare[1])))
            
            if gdat.verbtype > 1:
                print 'Element positions prior to calculating the pair list for the reverse move of a split'
                print 'lgal'
                print lgal * gdat.anglfact
                print 'bgal'
                print bgal * gdat.anglfact
    
    if gdatmodi.propmerg:
        
        # number of point sources to be modified
        gdatmodi.numbelemeval[0][0] = 3
        
        # determine the index of the primary element to be merged (in the full element list)
        gdatmodi.indxelemfullmergfrst = choice(arange(len(thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]])))

        ## first element index to be merged
        gdatmodi.mergindxelemfrst = thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmergfrst]

        # find the probability of merging this element with the others 
        probmerg = retr_probmerg(gdat, gdatmodi, 'this', gdatmodi.indxelemfullmergfrst)
        indxelemfulltemp = arange(len(thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]))
        gdatmodi.indxelemfullmergseco = choice(setdiff1d(indxelemfulltemp, array([gdatmodi.indxelemfullmergfrst])), p=probmerg)
        gdatmodi.indxelemfullmodi = sort(array([gdatmodi.indxelemfullmergfrst, gdatmodi.indxelemfullmergseco]))
        
        # parameters of the first element to be merged
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            ## first
            gdatmodi.compfrst[k] = thissampvarb[thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
        
        # determine indices of the modified elements in the sample vector
        ## first element
        # temp -- this would not work for multiple populations !
        gdatmodi.indxsampfrst = indxsampcompinit + numbcomp[gdatmodi.indxpoplmodi[0]] * gdatmodi.mergindxelemfrst
        gdatmodi.indxsampfrstfinl = gdatmodi.indxsampfrst + numbcomp[gdatmodi.indxpoplmodi[0]]
        gdatmodi.indxsamptran.append(arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl))

        ## second element index to be merged
        gdatmodi.mergindxelemseco = thisindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmergseco]
       
        ## second element
        gdatmodi.indxsampseco = indxsampcompinit + numbcomp[gdatmodi.indxpoplmodi[0]] * gdatmodi.mergindxelemseco
        gdatmodi.indxsampsecofinl = gdatmodi.indxsampseco + numbcomp[gdatmodi.indxpoplmodi[0]]
        gdatmodi.indxsamptran.append(arange(gdatmodi.indxsampseco, gdatmodi.indxsampsecofinl))
        
        # parameters of the elements to be merged
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            ## second
            gdatmodi.compseco[k] = thissampvarb[thisindxsampcomp[strgcomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[1]]]

        # indices of the sample vector elements to be modified
        gdatmodi.indxsampmodi = arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl)

        # indices of the element to be merged
        gdatmodi.indxelemmodi = [gdatmodi.mergindxelemfrst, gdatmodi.mergindxelemseco]

        # auxiliary parameters
        gdatmodi.auxipara[0] = gdatmodi.compseco[0] - gdatmodi.compfrst[0]
        gdatmodi.auxipara[1] = gdatmodi.compseco[1] - gdatmodi.compfrst[1]
        gdatmodi.auxipara[2] = gdatmodi.compfrst[2] / (gdatmodi.compfrst[2] + gdatmodi.compseco[2]) 
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            if k > 2:
                gdatmodi.auxipara[k] = gdatmodi.compseco[k]

        # merged element
        gdatmodi.comppare[2] = gdatmodi.compfrst[2] + gdatmodi.compseco[2]
        if gdatmodi.comppare[2] > getattr(gdat, 'maxm' + gdat.fittnamefeatampl[gdatmodi.indxpoplmodi[0]]):
            gdatmodi.thisaccpprio = False
            if gdat.verbtype > 1:
                print 'Proposal rejected due to falling outside the prior.'
        
        if False:
            for k in range(3):
                print 'k'
                print k
                print 'gdatmodi.compfrst[k]'
                print gdatmodi.compfrst[k]
                print 'gdatmodi.compseco[k]'
                print gdatmodi.compseco[k]
                print 'gdatmodi.comppare[k]'
                print gdatmodi.comppare[k]
                print
            if abs(gdatmodi.comppare[2] - 3e-8) / gdatmodi.comppare[2] > 0.001:
                raise Exception('')
            print 'gdatmodi.thisaccpprio'
            print gdatmodi.thisaccpprio
            print

        gdatmodi.comppare[0] = gdatmodi.compfrst[0] + (1. - gdatmodi.auxipara[2]) * (gdatmodi.compseco[0] - gdatmodi.compfrst[0])
        gdatmodi.comppare[1] = gdatmodi.compfrst[1] + (1. - gdatmodi.auxipara[2]) * (gdatmodi.compseco[1] - gdatmodi.compfrst[1])
        for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            if k < 2:
                gdatmodi.comppare[k] = gdatmodi.compfrst[k] + (1. - gdatmodi.auxipara[2]) * (gdatmodi.compseco[k] - gdatmodi.compfrst[k])
            elif k == 2:
                gdatmodi.comppare[k] = gdatmodi.compfrst[k] + gdatmodi.compseco[k]
            else:
                setattr(gdatmodi, strgcomp + 'pare', strgcomp + 'frst')

        gdatmodi.nextsamp[gdatmodi.indxsamptran[0][:3]] = cdfn_trap(gdat, gdatmodi, gdatmodi.comppare)[:3]
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0][3:]] = gdatmodi.thissamp[gdatmodi.indxsamptran[0][3:]]
        gdatmodi.nextsampvarb[gdatmodi.indxsamptran[0][:3]] = gdatmodi.comppare[:3]
        gdatmodi.nextsampvarb[gdatmodi.indxsamptran[0][3:]] = thissampvarb[gdatmodi.indxsamptran[0][3:]]

        # calculate the proposed list of pairs
        if gdat.verbtype > 1:
            print 'mergindxfrst: ', gdatmodi.mergindxelemfrst
            print 'gdatmodi.indxelemfullmergfrst: ', gdatmodi.indxelemfullmergfrst
            print 'mergindxseco: ', gdatmodi.mergindxelemseco
            print 'gdatmodi.indxelemfullmergseco: ', gdatmodi.indxelemfullmergseco
            print 'indxsampfrst: ', gdatmodi.indxsampfrst
            print 'gdatmodi.indxsampfrstfinl: ', gdatmodi.indxsampfrstfinl
            print 'indxsampseco: ', gdatmodi.indxsampseco
            print 'gdatmodi.indxsampsecofinl: ', gdatmodi.indxsampsecofinl
        
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
        gdatmodi.nextsamp[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] = \
                                                    gdatmodi.thissamp[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] + 1
        gdatmodi.nextsampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] = \
                                                    thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] + 1
    if gdatmodi.propdeth or gdatmodi.propmerg:
        nextsamp[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] = \
                                                                thissamp[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] - 1
        nextsampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] = \
                                                                thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]] - 1
    
    if (gdatmodi.propdeth or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        # remove the element from the occupied element list
        for a, indxelem in enumerate(gdatmodi.indxelemmodi):
            if a == 0 and gdatmodi.propdeth or a == 1 and gdatmodi.propmerg:
                nextindxelemfull[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]].remove(indxelem)
    
    if numbtrap > 0:
        if gdatmodi.propwith:
            if gdatmodi.propdist:
                gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampmodi]), \
                                        thisindxsampcomp[gdatmodi.thisstrgcomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]))
            gdatmodi.indxsampchec = gdatmodi.indxsampmodi
        else:
            gdatmodi.indxsampchec = []
            if (gdatmodi.propbrth or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
                gdatmodi.indxsampmodi = concatenate((indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0], None], gdatmodi.indxsamptran[0]))
            if gdatmodi.propdeth:
                gdatmodi.indxsampmodi = indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0], None]
            if gdatmodi.propsplt:
                gdatmodi.indxsampmodi = concatenate((indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0], None], \
                                                                                                                gdatmodi.indxsamptran[0], gdatmodi.indxsamptran[1]))
    else:
        gdatmodi.indxsampchec = gdat.indxfixpprop
        gdatmodi.indxsampmodi = gdat.indxfixpprop
    
    # reject the sample if proposal is outside the prior
    indxchecfail = where((nextsamp[gdatmodi.indxsampchec] < 0.) | (nextsamp[gdatmodi.indxsampchec] > 1.))[0]
    if indxchecfail.size > 0:
        if gdat.verbtype > 1:
            for k in range(20):
                print 'Proposal rejected due to proposal outside the prior during the common check'
            print 'indxchecfail'
            print indxchecfail
            print
        gdatmodi.thisaccpprio = False
    
    if False:
        print 'gdatmodi.thisaccpprio'
        print gdatmodi.thisaccpprio
        print
    
    if gdat.verbtype > 1:
        print 'gdatmodi.thisaccpprio'
        print gdatmodi.thisaccpprio
        print 'gdatmodi.indxsampchec'
        print gdatmodi.indxsampchec
        if numbtrap > 0:
            if gdatmodi.proptran:
                print 'gdatmodi.indxsamptran'
                for indxsamptran in gdatmodi.indxsamptran:
                    print indxsamptran
            print 'thisindxelemfull'
            print thisindxelemfull
            print 'nextindxelemfull'
            print nextindxelemfull
            if gdatmodi.propelem and not (gdatmodi.propmerg and not gdatmodi.thisaccpprio):
                print 'gdatmodi.indxelemmodi'
                print gdatmodi.indxelemmodi
                print 'gdatmodi.indxelemfullmodi'
                print gdatmodi.indxelemfullmodi

    if numbtrap > 0:
        if gdatmodi.propwith:
            nextindxsampcomp = thisindxsampcomp
        else:
            nextindxsampcomp = retr_indxsampcomp(gdat, nextindxelemfull, strgmodl)
        setattr(gdatobjt, strgpfixnext + 'indxsampcomp', nextindxsampcomp)
    
    # temp -- this is inefficient for propwithsing proposals
    if gdatmodi.thisaccpprio and (gdatmodi.propbrth):
        gdatmodi.nextsampvarb = retr_sampvarb(gdat, strgmodl, gdatmodi.nextsamp, gdatmodi.nextindxsampcomp)
   
    gdatmodi.dicteval = [[{}]]
    if gdatmodi.propelem:
        if boolelemspec[gdatmodi.indxpoplmodi[0]]:
		    gdatmodi.dicteval[0][0]['sind'] = []
		    gdatmodi.dicteval[0][0]['curv'] = []
		    gdatmodi.dicteval[0][0]['expc'] = []
		    gdatmodi.dicteval[0][0]['sind0001'] = []
		    gdatmodi.dicteval[0][0]['sind0002'] = []

        for k, namecomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
            if gdatmodi.propwith:
                thiscomp = thissampvarb[thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
                nextcomp = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
                if namecomp == namefeatampl[gdatmodi.indxpoplmodi[0]]:
                    gdatmodi.dicteval[0][0][namecomp] = array([-thiscomp, nextcomp])
                else:
                    gdatmodi.dicteval[0][0][namecomp] = array([thiscomp, nextcomp])
            elif gdatmodi.propbrth:
                comp = gdatmodi.nextsampvarb[gdatmodi.nextindxsampcomp[namecomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
                gdatmodi.dicteval[0][0][namecomp] = array([comp])
            elif gdatmodi.propdeth:
                comp = thissampvarb[thisindxsampcomp[namecomp][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]][gdatmodi.indxelemfullmodi[0]]]
                if namecomp == namefeatampl[gdatmodi.indxpoplmodi[0]]:
                    gdatmodi.dicteval[0][0][namecomp] = array([-comp])
                else:
                    gdatmodi.dicteval[0][0][namecomp] = array([comp])
            elif gdatmodi.propsplt:
                if namecomp == namefeatampl[gdatmodi.indxpoplmodi[0]]:
                    gdatmodi.dicteval[0][0][namecomp] = \
                                                                                           array([-gdatmodi.comppare[k], gdatmodi.compfrst[k], gdatmodi.compseco[k]])
                else:
                    gdatmodi.dicteval[0][0][namecomp] = array([gdatmodi.comppare[k], gdatmodi.compfrst[k], gdatmodi.compseco[k]])
            elif gdatmodi.propmerg and gdatmodi.thisaccpprio:
                if namecomp == namefeatampl[gdatmodi.indxpoplmodi[0]]:
                    gdatmodi.dicteval[0][0][namecomp] = array([-gdatmodi.compfrst[k], -gdatmodi.compseco[k], gdatmodi.comppare[k]])
                else:
                    gdatmodi.dicteval[0][0][namecomp] = array([gdatmodi.compfrst[k], gdatmodi.compseco[k], gdatmodi.comppare[k]])
        
        if boolelemspec[gdatmodi.indxpoplmodi[0]]:
            
            gdatmodi.dicteval[0][0]['spec'] = retr_spec(gdat, \
                                        gdatmodi.dicteval[0][0]['flux'], \
                                        sind=gdatmodi.dicteval[0][0]['sind'], \
                                        curv=gdatmodi.dicteval[0][0]['curv'], \
                                        expc=gdatmodi.dicteval[0][0]['expc'], \
                                        sind0001=gdatmodi.dicteval[0][0]['sind0001'], \
                                        sind0002=gdatmodi.dicteval[0][0]['sind0002'], \
                                        spectype=spectype[gdatmodi.indxpoplmodi[0]])
        if elemtype[gdatmodi.indxpoplmodi[0]] == 'lghtline':
            gdatmodi.dicteval[0][0]['spec'] = retr_spec(gdat, \
                                            gdatmodi.dicteval[0][0]['flux'], \
                                            elin=gdatmodi.dicteval[0][0]['elin'], edis=gdat.edis, spectype=spectype[gdatmodi.indxpoplmodi[0]])
            
    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        
        ## Jacobian
        if gdatmodi.propsplt:
            gdatmodi.thisljcbfact = log(gdatmodi.comppare[2])
        else:
            gdatmodi.thisljcbfact = -log(gdatmodi.comppare[2])
        
        thisnumbelem = thissampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]]
        nextnumbelem = gdatmodi.nextsampvarb[indxfixpnumbelem[gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]]
        
        ## combinatorial factor
        if gdatmodi.propmerg:
            probfrwdfrst = retr_probmerg(gdat, gdatmodi, 'this', gdatmodi.indxelemfullmodi[1], gdatmodi.indxelemfullmodi[0]) / thisnumbelem
            probfrwdseco = retr_probmerg(gdat, gdatmodi, 'this', gdatmodi.indxelemfullmodi[0], gdatmodi.indxelemfullmodi[1]) / thisnumbelem
            probreve = 1. / nextnumbelem
            gdatmodi.thislcomfact = log(probreve) - log(probfrwdfrst + probfrwdseco)
        else:
            probfrwd = 1. / thisnumbelem
            probrevefrst = retr_probmerg(gdat, gdatmodi, 'next', -1, gdatmodi.indxelemfullmodi[0]) / nextnumbelem
            probreveseco = retr_probmerg(gdat, gdatmodi, 'next', gdatmodi.indxelemfullmodi[0], -1) / nextnumbelem
            gdatmodi.thislcomfact = log(probrevefrst + probreveseco) - log(probfrwd)


def retr_probmerg(gdat, gdatmodi, strgstat, indxelemfullexcl, indxelemfulleval=None):
    
    lgal = getattr(gdatmodi, strgstat + 'sampvarb')[getattr(gdatmodi, strgstat + 'indxsampcomp')['lgal'][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]][indxelemfullexcl]
    bgal = getattr(gdatmodi, strgstat + 'sampvarb')[getattr(gdatmodi, strgstat + 'indxsampcomp')['bgal'][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]][indxelemfullexcl]
    
    lgalstat = getattr(gdatmodi, strgstat + 'sampvarb')[getattr(gdatmodi, strgstat + 'indxsampcomp')['lgal'][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]]
    bgalstat = getattr(gdatmodi, strgstat + 'sampvarb')[getattr(gdatmodi, strgstat + 'indxsampcomp')['bgal'][gdatmodi.indxpoplmodi[0]][gdatmodi.indxregimodi[0]]]
    
    angldist = retr_angldist(gdat, lgal, bgal, lgalstat, bgalstat)
    
    probmerg = exp((angldist - gdat.radispmr)**2)
    if indxelemfulleval == None:
        probmerg = concatenate((probmerg[:indxelemfullexcl], probmerg[indxelemfullexcl+1:]))
        probmerg /= sum(probmerg)
    else:
        if indxelemfullexcl == -1:
            probmerg = probmerg[:-1]
        else:
            probmerg = concatenate((probmerg[:indxelemfullexcl], probmerg[indxelemfullexcl+1:]))
        probmerg /= sum(probmerg)
        if indxelemfullexcl < indxelemfulleval and indxelemfulleval != -1:
            indxelemfulleval -= 1
        probmerg = probmerg[indxelemfulleval]

    return probmerg

    
def show_dbug(gdat, gdatmodi):
    print 'hey'
    print 'thissamp nextsamp diffsampvarb'
    for k in gdat.fittindxpara:
        print '%10.3g  %10.3g %10.3g' % (gdatmodi.thissamp[k], gdatmodi.nextsamp[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k])
    print

        
def retr_indxsamppnts(gdat, l, d, indxelem):

    indxsamppnts = gdat.fittindxsampcompinit + gdat.fittnumbtrapregipoplcuml[l][d] + indxelem[None, :] * gdat.fittnumbcomp[l] + gdat.fittindxcomp[l][:, None]
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
        fact = 2. * pi * trapz(psfnnorm * sin(gdat.binsangl[None, :, None]), gdat.binsangl, axis=1)[:, None, :]
        psfn /= fact

    return psfn


def retr_axis(gdat, strgvarb, minm=None, maxm=None, numb=None, bins=None, scal='self', invr=False, strginit=''):
    
    if bins == None:
        if scal == 'self' or scal == 'pois' or scal == 'gaus':
            binsscal = linspace(minm, maxm, numb + 1)
        if scal == 'logt':
            binsscal = linspace(log10(minm), log10(maxm), numb + 1)
        if scal == 'asnh':
            binsscal = linspace(arcsinh(minm), arcsinh(maxm), numb + 1)
        
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
            meanvarb = sinh(meanvarbscal)
            bins = sinh(binsscal)
    else:
        numb = bins.size - 1

    indx = arange(numb)
    delt = diff(bins) 
    limt = [amin(bins), amax(bins)] 

    setattr(gdat, strginit + 'limt' + strgvarb, limt)
    setattr(gdat, strginit + 'bins' + strgvarb, bins)
    setattr(gdat, strginit + 'mean' + strgvarb, meanvarb)
    setattr(gdat, strginit + 'delt' + strgvarb, delt)
    setattr(gdat, strginit + 'numb' + strgvarb, numb)
    setattr(gdat, strginit + 'indx' + strgvarb, indx)


def retr_unit(lgal, bgal):

    xdat = cos(bgal) * cos(lgal)
    ydat = -cos(bgal) * sin(lgal)
    zaxi = sin(bgal)

    return xdat, ydat, zaxi


def retr_psec(gdat, conv):

    # temp
    conv = conv.reshape((gdat.numbsidecart, gdat.numbsidecart))
    psec = (abs(fft.fft2(conv))**2)[:gdat.numbsidecart/2, :gdat.numbsidecart/2] * 1e-3
    psec = psec.flatten()

    return psec
   

def retr_psecodim(gdat, psec):
    
    psec = psec.reshape((gdat.numbsidecart / 2, gdat.numbsidecart / 2))
    psecodim = zeros(gdat.numbsidecart / 2)
    for k in gdat.indxmpolodim:
        indxmpol = where((gdat.meanmpol > gdat.binsmpolodim[k]) & (gdat.meanmpol < gdat.binsmpolodim[k+1]))
        psecodim[k] = mean(psec[indxmpol])
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
                #        indxelemcntr = indxtupl[indxstkstemp][2]
                #        gdatmodi.thissampvarb = gdat.listsampvarb[indxsamptotlcntr, :]
                #        gdatmodi.thisindxelemfull[].append()

                #    plot_genemaps(gdat, gdatmodi, 'this', 'cntpdata', indxenerplot=0, indxevttplot=0, cond=True)
                
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
            arry = stack(gdat.dictglob['liststkscond'][r][strgfeat], axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat] = zeros(([3] + list(arry.shape[1:])))
            gdat.dictglob['poststkscond'][r][strgfeat][0, ...] = median(arry, axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat][1, ...] = percentile(arry, 16., axis=0)
            gdat.dictglob['poststkscond'][r][strgfeat][2, ...] = percentile(arry, 84., axis=0)
            
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
    
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    # temp
    conv = abs(gradient(defl[:, :, 0], gdat.sizepixl, axis=0) + gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) / 2.
    conv = conv.flatten()
    return conv


def retr_invm(gdat, defl):
    
    # temp
    defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    invm = (1. - gradient(defl[:, :, 0], gdat.sizepixl, axis=0)) * (1. - gradient(defl[:, :, 1], gdat.sizepixl, axis=1)) - \
                                                                         gradient(defl[:, :, 0], gdat.sizepixl, axis=1) * gradient(defl[:, :, 1], gdat.sizepixl, axis=0)
    invm = invm.flatten()
    return invm


def setp_indxswepsave(gdat):

    gdat.indxswep = arange(gdat.numbswep)
    gdat.boolsave = zeros(gdat.numbswep, dtype=bool)
    gdat.indxswepsave = arange(gdat.numbburn, gdat.numbburn + gdat.numbsamp * gdat.factthin, gdat.factthin)
    gdat.boolsave[gdat.indxswepsave] = True
    gdat.indxsampsave = zeros(gdat.numbswep, dtype=int) - 1
    gdat.indxsampsave[gdat.indxswepsave] = arange(gdat.numbsamp)
    

def retr_cntspnts(gdat, lgal, bgal, spec, indxregipnts):
    
    indxpixlpnts = retr_indxpixl(gdat, bgal, lgal)
    cnts = zeros((gdat.numbener, lgal.size))
    for k in range(lgal.size):
        cnts[:, k] += spec[:, k] * gdat.expototl[indxregipnts[k]][:, indxpixlpnts[k]]
    if gdat.enerdiff:
        cnts *= gdat.deltener[:, None]
    cnts = sum(cnts, axis=0)

    return cnts


def retr_mdencrit(gdat, adissour, adishost, adishostsour):
    
    mdencrit = gdat.factnewtlght / 4. / pi * adissour / adishostsour / adishost
        
    return mdencrit


def retr_massfrombein(gdat, adissour, adishost, adishostsour):

    mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
    massfrombein = pi * adishost**2 * mdencrit

    return massfrombein


def retr_factmcutfromdefs(gdat, adissour, adishost, adishostsour, asca, acut):
    
    mdencrit = retr_mdencrit(gdat, adissour, adishost, adishostsour)
    
    fracacutasca = acut / asca
    
    factmcutfromdefs = pi * adishost**2 * mdencrit * asca * retr_mcutfrommscl(fracacutasca)

    return factmcutfromdefs


def retr_mcut(gdat, defs, asca, acut):
    
    mscl = defs * pi * gdat.adishost**2 * gdat.mdencrit * asca
    fracacutasca = acut / asca
    mcut = mscl * retr_mcutfrommscl(fracacutasca)
    
    return mcut


def retr_mcutfrommscl(fracacutasca):
    
    mcut = fracacutasca**2 / (fracacutasca**2 + 1.)**2 * ((fracacutasca**2 - 1.) * log(fracacutasca) + fracacutasca * pi - (fracacutasca**2 + 1.))

    return mcut


def retr_negalogt(varb):
    
    negalogt = sign(varb) * log10(fabs(varb))
    
    return negalogt


def setpprem(gdat):

    # set up the indices of the models
    for strgmodl in gdat.liststrgmodl:
        retr_indxsamp(gdat, strgmodl=strgmodl, init=True)
    
    gdat.liststrglimt = ['minm', 'maxm']
   
    gdat.indxregi = arange(gdat.numbregi)

    # temp
    gdat.edis = 1. / 2.35

    gdat.lablcurv = r'\kappa'
    gdat.lablexpc = r'E_{c}'
    
    gdat.scalcurvplot = 'self'
    gdat.scalexpcplot = 'self'
    gdat.factcurvplot = 1.
    gdat.factexpcplot = 1.
    
    if gdat.exprtype == 'hubb':
        gdat.numbgrid = 3
    else:
        gdat.numbgrid = 1
    gdat.indxgrid = arange(gdat.numbgrid)

    # axes
    gdat.minmlgaldata = -gdat.maxmgangdata
    gdat.maxmlgaldata = gdat.maxmgangdata
    gdat.minmbgaldata = -gdat.maxmgangdata
    gdat.maxmbgaldata = gdat.maxmgangdata
    
    if gdat.pixltype == 'heal' and gdat.forccart:
        raise Exception('Cartesian forcing can only used with cart pixltype')

    # input data
    gdat.pathinpt = gdat.pathdata + 'inpt/'
    if gdat.datatype == 'inpt':
        path = gdat.pathinpt + gdat.strgexprsbrt
        gdat.sbrtdata = pf.getdata(path)
        
        # temp
        if (gdat.pixltype == 'heal' or gdat.pixltype == 'cart' and gdat.forccart) and gdat.sbrtdata.ndim == 3 or gdat.pixltype == 'cart' and gdat.sbrtdata.ndim == 4:
            print 'Input data incompatible with PCAT %s. Converting...' % gdat.strgvers
            gdat.sbrtdata = gdat.sbrtdata[None, :, :, :]
        
        # temp
        if gdat.exprtype == 'ferm' and gdat.sbrtdata.shape[3] == 4 and (gdat.anlytype.startswith('rec7') or gdat.anlytype.startswith('manu')):
            print 'Input data incompatible with new Fermi-LAT event type binning. Converting...'
            gdat.sbrtdata = gdat.sbrtdata[:, :, :, 2:4]

        if gdat.pixltype == 'heal' or gdat.pixltype == 'cart' and gdat.forccart:
            if gdat.sbrtdata.ndim != 4:
                raise Exception('exprsbrtdata should be a 4D numpy array if pixelization is HealPix.')
        else:
            if gdat.sbrtdata.ndim != 5:
                raise Exception('exprsbrtdata should be a 5D numpy array if pixelization is Cartesian.')
        
        if gdat.pixltype == 'cart':
            if gdat.forccart:
                gdat.sbrtdatatemp = empty((gdat.sbrtdata.shape[0], gdat.sbrtdata.shape[1], gdat.numbsidecart, gdat.numbsidecart, gdat.sbrtdata.shape[3]))
                for d in arange(gdat.sbrtdata.shape[0]):
                    for i in arange(gdat.sbrtdata.shape[1]):
                        for m in arange(gdat.sbrtdata.shape[3]):
                            gdat.sbrtdatatemp[d, i, :, :, m] = tdpy.util.retr_cart(gdat.sbrtdata[d, i, :, m], numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).T
                gdat.sbrtdata = gdat.sbrtdatatemp
            else:
                gdat.numbsidecart = gdat.sbrtdata.shape[2]
            gdat.sbrtdata = gdat.sbrtdata.reshape((gdat.sbrtdata.shape[0], gdat.sbrtdata.shape[1], gdat.numbsidecart**2, gdat.sbrtdata.shape[4]))
        
        gdat.numbenerfull = gdat.sbrtdata.shape[1]
        gdat.numbpixlfull = gdat.sbrtdata.shape[2]
        gdat.numbevttfull = gdat.sbrtdata.shape[3]

        if gdat.pixltype == 'heal':
            # temp
            gdat.numbsidecart = 100
            gdat.numbsideheal = int(sqrt(gdat.numbpixlfull / 12))
    
    gdat.minmcurv = -1.
    gdat.maxmcurv = 1.
    gdat.minmexpc = 0.1
    gdat.maxmexpc = 10.
        
    # pixelization
    if gdat.pixltype == 'cart':
        gdat.apix = (2. * gdat.maxmgangdata / gdat.numbsidecart)**2
    if gdat.pixltype == 'heal':
        lgalheal, bgalheal, gdat.numbpixlfull, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)
    gdat.sizepixl = sqrt(gdat.apix)
    
    if gdat.datatype == 'mock':
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2
        if gdat.pixltype == 'heal':
            gdat.numbpixlfull = 12 * gdat.numbsideheal**2
    

def setpinit(gdat, boolinitsetp=False):

    if False and gdat.elemtype == 'lens' and gdat.strgproc == 'fink2.rc.fas.harvard.edu':
        cliblens = ctypes.CDLL(os.environ["PCAT_PATH"] + '/cliblens.so')
        cliblens.retr_deflsubh()
    
    # set up the indices of the fitting model
    retr_indxsamp(gdat)
    
    gdat.listnamevarbstat = ['samp', 'sampvarb', 'indxelemfull', 'indxsampcomp', 'lliktotl', 'llik', 'lpritotl', 'lpri']
    for name in gdat.fittlistnamediff:
        #if not (name.startswith('back') and gdat.fittunifback[int(name[4:])]):
        gdat.listnamevarbstat += ['sbrt' + name + 'conv']
    if True in gdat.fittboolelemfore:
        gdat.listnamevarbstat += ['sbrtpnts']
    if gdat.fittlensmodltype != 'none':
        gdat.listnamevarbstat += ['deflhost']
    if gdat.fitthostemistype != 'none':
        gdat.listnamevarbstat += ['sbrthost']
    if gdat.fittpsfnevaltype == 'conv' or gdat.fittpsfnevaltype == 'full':
        gdat.listnamevarbstat += ['psfnconv']
    if 'lens' in gdat.fittelemtype:
        gdat.listnamevarbstat += ['deflelem']
    
    # paths
    ## data
    gdat.pathpixlcnvt = gdat.pathdata + 'pixlcnvt/'
    gdat.pathprox = gdat.pathdata + 'prox/'
    ## plot
    gdat.pathplot = gdat.pathimag + gdat.rtag + '/'
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
            
        for namethrd in ['hist', 'trac', 'join', 'cova']:
            setattr(gdat, 'path' + name + 'finlvarbscal' + namethrd, path + 'finl/varbscal/' + namethrd + '/')
            
        for nameseco in ['histodim', 'histtdim', 'assc', 'scattdim', 'cmpl', 'fdis']:
            for namethrd in ['init', 'fram', 'finl', 'anim']:
                if namethrd == 'init':
                    if nameseco == 'assc' or nameseco == 'histtdim' or nameseco == 'fdis' or nameseco == 'cmpl':
                        continue
                    setattr(gdat, 'path' + namethrd + nameseco, gdat.pathplot + 'init/' + nameseco + '/')
                else:
                    setattr(gdat, 'path' + name + namethrd + nameseco, path + namethrd + '/' + nameseco + '/')
            
    gdat.pathinfo = gdat.pathplot + 'info/'
    
    ## make the directories 
    for attr, valu in gdat.__dict__.iteritems():
        if attr.startswith('path'):
            os.system('mkdir -p %s' % valu)
 
    ## names of the variables for which cumulative posteriors will be plotted
    if 'lens' in gdat.fittelemtype:
        gdat.listnamevarbcpos = ['convelem']
    else:
        gdat.listnamevarbcpos = []

    gdat.ascaglob = 0.05 / gdat.anglfact
    gdat.acutglob = 1. / gdat.anglfact
    gdat.cutfdefs = 3e-3 / gdat.anglfact

    # plotting
    gdat.legdsampdist = 'Posterior'
    gdat.legdsamp = 'Sample'
    gdat.legdmlik = 'Maximum likelihood'
    gdat.legdmedi = 'Median'
    gdat.legdstdv = 'Std. dev.'
    
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
    for d in gdat.indxregi:
        strgregi = 'reg%d' % d
        setattr(gdat, 'lablfracsubh' + strgregi, r'$f_{\rm{sub}}$')
        setattr(gdat, 'lablfracsubhdeltbein' + strgregi, r'$f_{\rm{sub,E,%d}}$' % d)
        setattr(gdat, 'lablfracsubhintgbein' + strgregi, r'$f_{\rm{sub,E,%d<}}$' % d)
        setattr(gdat, 'lablmasssubhdeltbein' + strgregi, r'$M_{\rm{sub,E,%d}}$' % d)
        setattr(gdat, 'lablmasssubhintgbein' + strgregi, r'$M_{\rm{sub,E,%d<}}$' % d)
        setattr(gdat, 'lablmasshostdeltbein' + strgregi, r'$M_{\rm{hst,E,%d}}$' % d)
        setattr(gdat, 'lablmasshostintgbein' + strgregi, r'$M_{\rm{hst,E,%d<}}$' % d)
        for namevarb in ['fracsubh', 'masssubh', 'masshost']:
            for namecalc in ['delt', 'intg']:
                for nameeval in ['', 'bein']:
                    setattr(gdat, 'scal' + namevarb + namecalc + nameeval + strgregi, 'logt')
        setattr(gdat, 'lablmassintg' + strgregi, r'<M>_{<r}')
        setattr(gdat, 'lablmassintgunit' + strgregi, r'$M_{\odot}$')
        setattr(gdat, 'lablmassdelt' + strgregi, r'<M>_r')
        setattr(gdat, 'lablmassdeltunit' + strgregi, r'$M_{\odot}$')
        setattr(gdat, 'lablfracsubhintg' + strgregi, r'<f>_{\rm{<r,sub}}')
        setattr(gdat, 'lablfracsubhdelt' + strgregi, r'<f>_{\rm{r,sub}}')
    
    gdat.minmreds = 0.
    gdat.maxmreds = 8.
    
    gdat.minmmagt = 20.
    gdat.maxmmagt = 30.

    gdat.scalmaxmnumbelem = 'logt'
    gdat.scalmedilliktotl = 'logt'

    gdat.lablener = 'E'
    #gdat.lablenertotl = '$%s$ [%s]' % (gdat.lablener, gdat.strgenerunit)
    
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

    gdat.labldefsunit = u'$^{\prime\prime}$'
    gdat.lablprat = 'cm$^{-2}$ s$^{-1}$'
    for nameenerscaltype in ['ene0', 'ene1', 'ene2', 'ene3']:
        
        for labltemptemp in ['flux', 'sbrt']:
            labltemp = getattr(gdat, 'labl' + labltemptemp)

            # define the label
            if nameenerscaltype == 'ene0':
                strgenerscal = '%s' % labltemp
            if nameenerscaltype == 'ene1':
                strgenerscal = 'E%s' % labltemp
            if nameenerscaltype == 'ene2':
                strgenerscal = 'E^2%s' % labltemp
            if nameenerscaltype == 'ene3':
                strgenerscal = '%s' % labltemp
            labl = '%s' % strgenerscal
            setattr(gdat, 'labl' + labltemptemp + nameenerscaltype, labl)

            for nameenerunit in ['gevv', 'ergs', 'kevv']:
                
                strgenerunit = getattr(gdat, 'strgener' + nameenerunit)

                if nameenerscaltype == 'ene0':
                    strgenerscalunit = '%s$^{-1}$' % strgenerunit
                if nameenerscaltype == 'ene1':
                    strgenerscalunit = '' 
                if nameenerscaltype == 'ene2':
                    strgenerscalunit = '%s' % strgenerunit
                if nameenerscaltype == 'ene3':
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
        gdat.lablfluxunit = getattr(gdat, 'lablfluxene0' + gdat.nameenerunit + 'unit')
        gdat.lablsbrtunit = getattr(gdat, 'lablsbrtene0' + gdat.nameenerunit + 'sterunit')

    gdat.lablexpo = r'$\epsilon$'
    gdat.lablexpounit = 'cm$^2$ s'
    
    gdat.lablprvl = '$p$'
    
    gdat.lablreds = 'z'
    gdat.lablmagt = 'm_R'

    gdat.lablsind = 's'
    gdat.lablsind0001 = 's_1'
    gdat.lablsind0002 = 's_2'
    gdat.lablexpcunit = gdat.strgenerunit
    
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
   
    for q in gdat.indxrefr:
        for l in gdat.fittindxpopl:
            for d in gdat.fittindxregipopl[l]:
                setattr(gdat, 'minmfdisref%dpop%dreg%d' % (q, l, d), 0.)
                setattr(gdat, 'maxmfdisref%dpop%dreg%d' % (q, l, d), 1.)
                setattr(gdat, 'corrfdisref%dpop%dreg%d' % (q, l, d), 1.)
                setattr(gdat, 'factfdisref%dpop%dreg%dplot' % (q, l, d), 1.)
                setattr(gdat, 'scalfdisref%dpop%dreg%d' % (q, l, d), 'self')
                setattr(gdat, 'lablfdisref%dpop%dreg%d' % (q, l, d), '$f_{%d%d%d}$' % (q, l, d))
                
                setattr(gdat, 'minmcmplref%dpop%dreg%d' % (q, l, d), 0.)
                setattr(gdat, 'maxmcmplref%dpop%dreg%d' % (q, l, d), 1.)
                setattr(gdat, 'corrcmplref%dpop%dreg%d' % (q, l, d), 1.)
                setattr(gdat, 'factcmplref%dpop%dreg%dplot' % (q, l, d), 1.)
                setattr(gdat, 'scalcmplref%dpop%dreg%d' % (q, l, d), 'self')
                setattr(gdat, 'lablcmplref%dpop%dreg%d' % (q, l, d), '$c_{%d%d%d}$' % (q, l, d))
    
    dicttemp = deepcopy(gdat.__dict__)
    for name, valu in dicttemp.iteritems():
        if name.startswith('labl') and not name.endswith('unit'):
            name = name[4:]
            labl = getattr(gdat, 'labl' + name)
            try:
                lablunit = getattr(gdat, 'labl' + name + 'unit')
                if lablunit == '' or lablunit == None:
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
    
    # temp
    gdat.maxmgang = max(gdat.fittmaxmgang, gdat.truemaxmgang)
    
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
        gdat.minmcnts = 1e1
        gdat.maxmcnts = 1e5
    if gdat.exprtype == 'chan':
        gdat.minmcnts = 1.
        gdat.maxmcnts = 1e3
    if gdat.exprtype == 'hubb':
        gdat.minmcnts = 1.
        gdat.maxmcnts = 1e3

    gdat.minmspecplot = gdat.minmspec
    gdat.maxmspecplot = gdat.maxmspec
    
    gdat.minmdeltllik = 1.
    gdat.maxmdeltllik = 1e3
    gdat.minmdiss = 0.
    gdat.maxmdiss = gdat.maxmgang * sqrt(2.)
    
    gdat.minmrele = 0.01
    gdat.maxmrele = 1.

    gdat.minmreln = 0.3
    gdat.maxmreln = 10.

    gdat.minmreld = 0.1
    gdat.maxmreld = 10.

    gdat.minmrelc = 0.01
    gdat.maxmrelc = 1.

    gdat.minmmcut = 3e7
    gdat.maxmmcut = 2e9
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

    gdat.minmmasshost = 1e10
    gdat.maxmmasshost = 1e13
    gdat.scalmasshost = 'self'
    
    gdat.minmmasssubh = 1e8
    gdat.maxmmasssubh = 1e10
    gdat.scalmasssubh = 'self'

    if gdat.datatype == 'inpt':
        for l in gdat.fittindxpopl:
            for strgpdfn in gdat.fittliststrgpdfnprio[l]:
                if strgpdfn.startswith('gaum') and gdat.fittlgalprio == None and gdat.fittbgalprio == None:
                    raise Exception('If spatdisttype is "gaus", spatial coordinates of the prior catalog should be provided via lgalprio and bgalprio.')
    
    if gdat.propcomp == None:
        gdat.propcomp = gdat.fittnumbtrap > 0

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
                if strgfeat.startswith('defs') or strgfeat == 'gang' or strgfeat == 'lgal' or strgfeat == 'bgal' or \
                                                                                                            strgfeat == 'diss' or strgfeat == 'asca' or strgfeat == 'acut':
                    setattr(gdat, 'fact' + strgfeat + 'plot', gdat.anglfact)
                else:
                    setattr(gdat, 'fact' + strgfeat + 'plot', 1.)
                
                if strgfeat.startswith(gdat.namefeatsign) or strgfeat.startswith(namefeatampl[l]) or strgfeat == 'expo' or \
                                     strgfeat == 'cnts' or strgfeat.startswith('per') or \
                                     strgfeat == 'relc' or strgfeat == 'rele' or strgfeat == 'reln' or strgfeat == 'reld' or strgfeat.startswith('mcut') or strgfeat == 'deltllik':
                    setattr(gdat, 'scal' + strgfeat + 'plot', 'logt')
                else:
                    setattr(gdat, 'scal' + strgfeat + 'plot', 'self')
    
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
    
    if gdat.datatype != 'mock':
        if True in gdat.fittboolelemspec and gdat.exprtype == 'ferm' and gdat.maxmgangdata == 20. / gdat.anglfact:
            path = gdat.pathinpt + 'sbrt0018.png'
            gdat.sbrt0018 = sp.ndimage.imread(path, flatten=True)
            gdat.sbrt0018 -= amin(gdat.sbrt0018)
            gdat.sbrt0018 /= amax(gdat.sbrt0018)
            binslgaltemp = linspace(-gdat.fittmaxmgang, gdat.fittmaxmgang, gdat.sbrt0018.shape[1])
            binsbgaltemp = linspace(-gdat.fittmaxmgang, gdat.fittmaxmgang, gdat.sbrt0018.shape[0])
            gdat.sbrt0018objt = sp.interpolate.RectBivariateSpline(binsbgaltemp, binslgaltemp, gdat.sbrt0018)

    # list of scalar variable names
    gdat.listnamevarbscal = list(gdat.fittnamefixp)
    booltemp = False
    for strgmodl in gdat.liststrgmodl:
        lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
        if lensmodltype != 'none':
            booltemp = True
    if booltemp:
        for d in gdat.indxregi:
            strgregi = 'reg%d' % d
            for strgtemp in ['delt', 'intg']:
                gdat.listnamevarbscal += ['masshost' + strgtemp + 'bein' + strgregi]
                setattr(gdat, 'minmmasshost' + strgtemp + 'bein' + strgregi, gdat.minmmasshost)
                setattr(gdat, 'maxmmasshost' + strgtemp + 'bein' + strgregi, gdat.maxmmasshost)
            if gdat.fittnumbtrap > 0:
                if 'lens' in gdat.fittelemtype:
                    for strgtemp in ['delt', 'intg']:
                        gdat.listnamevarbscal += ['masssubh' + strgtemp + 'bein' + strgregi, 'fracsubh' + strgtemp + 'bein' + strgregi] 
                        setattr(gdat, 'minmmasssubh' + strgtemp + 'bein' + strgregi, gdat.minmmasssubh)
                        setattr(gdat, 'maxmmasssubh' + strgtemp + 'bein' + strgregi, gdat.maxmmasssubh)
                        setattr(gdat, 'minmfracsubh' + strgtemp + 'bein' + strgregi, gdat.minmfracsubh)
                        setattr(gdat, 'maxmfracsubh' + strgtemp + 'bein' + strgregi, gdat.maxmfracsubh)
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
    
    # log-prior register
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    numb = 0
    for l in gdat.fittindxpopl:
        numb += len(gdat.fittliststrgfeat[l])
    if gdat.fittmaxmnumbelemtotl > 0:
        if gdat.datatype == 'mock':
            gdat.numblpri = 3 + max(gdat.fittmaxmnumbcomp * gdat.fittnumbpopl, gdat.truemaxmnumbcomp * gdat.truemaxmnumbcomp) * gdat.numbregi
        else:
            gdat.numblpri = 3 + gdat.fittmaxmnumbcomp * gdat.fittnumbpopl * gdat.numbregi
    else:
        gdat.numblpri = 0
    if gdat.penalpridiff:
        gdat.numblpri += 1

    # size of the auxiliary variable propobability density vector
    gdat.numblpau = gdat.fittmaxmnumbcomp
    
    # process index
    gdat.indxproc = arange(gdat.numbproc)

    fileh5py = h5py.File(gdat.pathdata + 'inpt/adis.h5','r')
    reds = fileh5py['reds'][()]
    adis = fileh5py['adis'][()]
    adistdim = fileh5py['adistdim'][()]
    gdat.adisobjt = interp1d_pick(reds, adis)
    fileh5py.close()

    if 'lens' in gdat.commelemtype:
        gdat.redshost = 0.2
        gdat.redssour = 1.
        gdat.adishost = gdat.adisobjt(gdat.redshost) * 1e3 # [kpc]
        gdat.adissour = gdat.adisobjt(gdat.redssour) * 1e3 # [kpc]
        gdat.adishostsour = gdat.adissour - (1. + gdat.redshost) / (1. + gdat.redssour) * gdat.adishost
        gdat.adisfact = gdat.adishost * gdat.adissour / gdat.adishostsour
        gdat.factnewtlght = 2.09e16 # Msun / kpc
        gdat.massfrombein = retr_massfrombein(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        gdat.mdencrit = retr_mdencrit(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        
        retr_axis(gdat, 'mcut', gdat.minmmcut, gdat.maxmmcut, gdat.numbbinsplot)
        retr_axis(gdat, 'bein', gdat.minmbein, gdat.maxmbein, gdat.numbbinsplot)

    # angular deviation
    gdat.numbanglhalf = 10
    gdat.indxanglhalf = arange(gdat.numbanglhalf)
    retr_axis(gdat, 'anglhalf', 0., gdat.maxmgangdata, gdat.numbanglhalf)
    gdat.numbanglfull = 1000
    gdat.indxanglfull = arange(gdat.numbanglfull)
    retr_axis(gdat, 'anglfull', 0., 3. * gdat.maxmgang, gdat.numbanglfull)
    
    # temp
    #gdat.binsanglcosi = sort(cos(gdat.binsangl))
    
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
    if gdat.allwrefr:
        if gdat.datatype == 'mock':
            if gdat.truelegdpopl != None:
                # user-defined true population legend
                gdat.legdrefr = gdat.truelegdpopl
            else:
                # default true population legend
                gdat.legdrefr = []
                for l in gdat.trueindxpopl:
                    if gdat.truenumbpopl == 1:
                        gdat.legdrefr.append('True')
                    else:
                        gdat.legdrefr.append('True Pop. %d' % l)
        gdat.legdrefrmiss = []
        gdat.legdrefrhits = []
        for q in gdat.indxrefr:
            gdat.legdrefrmiss.append(gdat.legdrefr[q] + ' miss')
            gdat.legdrefrhits.append(gdat.legdrefr[q] + ' hit')
        if gdat.datatype == 'mock':
            gdat.legdrefrhost = gdat.legdrefr[0] + ' host'
            gdat.legdrefrsour = gdat.legdrefr[0] + ' sour'
   
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

    gdat.listnamechro = ['totl', 'type', 'prop', 'diag', 'save', 'plot', 'proc', 'lpri', 'llik', 'sbrtmodl', 'sbrtdiffconv']
    gdat.listlegdchro = ['Total', 'Type', 'Proposal', 'Diagnostics', 'Save', 'Plot', 'Process', 'Prior', 'Posterior', 'Total emission', 'Diffuse Conv.']
    if gdat.fittnumbtrap > 0:
        gdat.listnamechro += ['evalelem', 'kernelem']
        gdat.listlegdchro += ['El. Eval. Arrays', 'Kern. Ev.']
    if 'lens' in gdat.fittelemtype:
        gdat.listnamechro += ['deflzero', 'deflhost', 'deflextr', 'sbrtlens', 'sbrthost']
        gdat.listlegdchro += ['Array initialization', 'Host Deflection', 'External deflection', 'Lensed emission', 'Host emission']
    
    booltemp = False
    for strgmodl in gdat.liststrgmodl:
        booltemp = booltemp or getattr(gdat, strgmodl + 'psfnevaltype')
    if booltemp or gdat.fittpsfnevaltype == 'full' or gdat.truepsfnevaltype == 'full':
        gdat.listnamechro += ['psfnconv']
        gdat.listlegdchro += ['Img for PSF Conv.']
    
    gdat.listnamechro += ['expo', 'llikcalc']
    gdat.listlegdchro += ['Exposure', 'Log-likelihood']
    
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

    ## longitude
    gdat.numblgalpntsprob = gdat.numbsidepntsprob
    gdat.numbbgalpntsprob = gdat.numbsidepntsprob
    gdat.binslgalpntsprob = linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.binsbgalpntsprob = linspace(-gdat.maxmgangdata, gdat.maxmgangdata, gdat.numbsidepntsprob + 1)
    gdat.indxlgalpntsprob = arange(gdat.numblgalpntsprob)
    gdat.indxbgalpntsprob = arange(gdat.numbbgalpntsprob)

    retr_axis(gdat, 'defl', -gdat.maxmgang, gdat.maxmgang, 40)
    retr_axis(gdat, 'deflelem', -gdat.maxmgang * 1e-2, gdat.maxmgang * 1e-2, 40)

    # lensing problem setup
    ## number of deflection components to plot

    gdat.binslgalcartmesh, gdat.binsbgalcartmesh = meshgrid(gdat.binslgalcart, gdat.binsbgalcart)
    gdat.meanlgalcartmesh, gdat.meanbgalcartmesh = meshgrid(gdat.meanlgalcart, gdat.meanbgalcart)
    if gdat.pixltype == 'cart':
        gdat.sizepixl = sqrt(gdat.apix)
        gdat.indxpixlrofi = arange(gdat.numbsidecart**2)
        gdat.indxsidecart = arange(gdat.numbsidecart)
        gdat.indxsidemesh = meshgrid(gdat.indxsidecart, gdat.indxsidecart, indexing='ij')
        gdat.lgalgrid = gdat.meanlgalcart[gdat.indxsidemesh[0].flatten()]
        gdat.bgalgrid = gdat.meanbgalcart[gdat.indxsidemesh[1].flatten()]
        gdat.shapcart = (gdat.numbsidecart, gdat.numbsidecart)
        gdat.lgalgridcart = gdat.lgalgrid.reshape(gdat.shapcart)
        gdat.bgalgridcart = gdat.bgalgrid.reshape(gdat.shapcart)
        gdat.indxpent = meshgrid(gdat.indxregi, gdat.indxener, gdat.indxsidecart, gdat.indxsidecart, gdat.indxevtt, indexing='ij')
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
    if 'lens' in gdat.commelemtype:
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
        if 'lens' in gdat.commelemtype:
            gdat.minmwvecodim = gdat.minmmpolodim / gdat.adis
            gdat.maxmwvecodim = gdat.maxmmpolodim / gdat.adis
            gdat.minmwlenodim = gdat.minmanglodim * gdat.adis
            gdat.maxmwlenodim = gdat.maxmanglodim * gdat.adis
            retr_axis(gdat, 'wvecodim', gdat.minmwvecodim, gdat.maxmwvecodim, gdat.numbsidecart / 2)
            retr_axis(gdat, 'wlenodim', gdat.minmwlenodim, gdat.maxmwlenodim, gdat.numbsidecart, invr=True)
            gdat.meanwveclgal, gdat.meanwvecbgal = meshgrid(gdat.meanwvecodim, gdat.meanwvecodim, indexing='ij')
            gdat.meanwvec = sqrt(gdat.meanwveclgal**2 + gdat.meanwvecbgal**2)
        gdat.meanmpollgal, gdat.meanmpolbgal = meshgrid(gdat.meanmpolodim, gdat.meanmpolodim, indexing='ij')
        gdat.meanmpol = sqrt(gdat.meanmpollgal**2 + gdat.meanmpolbgal**2)

    # element parameter vector indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompflux = 2
    gdat.indxcompsind = 3
    gdat.indxcompcurv = 4
    gdat.indxcompexpc = 4

    # exposure
    if gdat.correxpo:
        if isinstance(gdat.strgexpo, float):
            if gdat.datatype == 'mock':
                if gdat.pixltype == 'heal':
                    gdat.expo = gdat.strgexpo * ones((gdat.numbregi, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                if gdat.pixltype == 'cart':
                    gdat.expo = gdat.strgexpo * ones((gdat.numbregi, gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
            if gdat.datatype == 'inpt':
                gdat.expo = gdat.strgexpo * ones((gdat.numbregi, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
        else: 
            path = gdat.pathinpt + gdat.strgexpo
            gdat.expo = pf.getdata(path)
            
            # temp
            if (gdat.pixltype == 'heal' or gdat.pixltype == 'cart' and gdat.forccart) and gdat.expo.ndim == 3 or gdat.pixltype == 'cart' and gdat.expo.ndim == 4:
                print 'Exposure map incompatible with PCAT %s. Converting...' % gdat.strgvers
                gdat.expo = gdat.expo[None, :, :, :]

            # temp
            if gdat.exprtype == 'ferm' and gdat.expo.shape[3] == 4 and (gdat.anlytype.startswith('rec7') or gdat.anlytype.startswith('manu')):
                print 'Input exposure incompatible with new Fermi-LAT event type binning. Converting...'
                gdat.expo = gdat.expo[:, :, :, 2:4]

            if amin(gdat.expo) == amax(gdat.expo):
                raise Exception('Bad input exposure map.')
            if gdat.pixltype == 'cart':
   
                if gdat.forccart:
                    gdat.expotemp = empty((gdat.numbregi, gdat.numbenerfull, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevttfull))
                    for d in gdat.indxregi:
                        for i in gdat.indxenerfull:
                            for m in gdat.indxevttfull:
                                gdat.expotemp[d, i, :, :, m] = tdpy.util.retr_cart(gdat.expo[d, i, :, m], numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                                                             minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                                             minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).T
                    gdat.expo = gdat.expotemp
                
                gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], gdat.expo.shape[1], -1, gdat.expo.shape[-1]))
    
    if gdat.killexpo:
        gdat.expo *= 1e-10
    if gdat.highexpo:
        gdat.expo *= 1e10

    # only include desired energy and PSF class bins
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    
    ## exposure
    if gdat.correxpo:
        # temp -- for some reason lists of arrays require manual processing
        gdat.expo = rplc_list(gdat.expo, gdat.indxcubeincl)
        if gdat.datatype == 'inpt':
            gdat.sbrtdata = rplc_list(gdat.sbrtdata, gdat.indxcubeincl)

    ## backgrounds
    for strgmodl in gdat.liststrgmodl:
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        indxback = getattr(gdat, strgmodl + 'indxback')
        sbrtbacknormincl = [[[] for d in gdat.indxregi] for c in indxback]
        for c in indxback:
            for d in gdat.indxregi:
                sbrtbacknormincl[c][d] = sbrtbacknorm[c][d][gdat.indxcubeincl]
        setattr(gdat, strgmodl + 'sbrtbacknorm', sbrtbacknormincl)
    
    # obtain cartesian versions of the maps
    #if gdat.pixltype == 'cart':
    #    gdat.expocart = gdat.expo.reshape((gdat.numbregi, gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
    #    for strgmodl in gdat.liststrgmodl:
    #        sbrtbacknormcart = []
    #        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
    #        for c in getattr(gdat, strgmodl + 'indxback'):
    #            sbrtbacknormcart.append(sbrtbacknorm[c].reshape((gdat.numbregi, gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt)))
    #        setattr(gdat, strgmodl + 'sbrtbacknormcart', sbrtbacknormcart)
    
    # mask the exposure map
    if gdat.listmask != None:
        for mask in gdat.listmask:
            if mask[0] == 'sqre':
                indxpixlmask = where((gdat.lgalgrid > mask[1]) & (gdat.lgalgrid < mask[2]) & (gdat.bgalgrid > mask[3]) & (gdat.bgalgrid < mask[4]))[0]
            if mask[0] == 'circ':
                indxpixlmask = where(sqrt((gdat.lgalgrid - mask[1])**2 + (gdat.bgalgrid - mask[2])**2) < mask[3])[0]
            for d in gdat.indxregi:
                gdat.expo[d][:, indxpixlmask, :] = 0.

    # plotting
    ## ROI
    gdat.exttrofi = array([gdat.minmlgaldata, gdat.maxmlgaldata, gdat.minmbgaldata, gdat.maxmbgaldata])
    gdat.exttrofi *= gdat.anglfact 
    gdat.frambndrdata = gdat.maxmgangdata * gdat.anglfact
    gdat.frambndrmodl = gdat.maxmgang * gdat.anglfact

    ## marker size
    gdat.minmmrkrsize = 100
    gdat.maxmmrkrsize = 500
    ## marker line width
    gdat.mrkrlinewdth = 3
    ## marker opacity
    gdat.alphmrkr = 0.5
    gdat.alphbndr = 0.5
    gdat.alphelem = 0.6
    gdat.alphmaps = 1.
    
    # number of colorbar ticks in the maps
    gdat.numbtickcbar = 11
    
    ## color bars
    gdat.minmlpdfspatpriointp = log(1. / 2. / gdat.maxmgang) - 10.
    gdat.maxmlpdfspatpriointp = log(1. / 2. / gdat.maxmgang) + 10.
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
    gdat.cmapcntpresi = 'RdBu'
    
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
   
    if amax(gdat.expo) <= 0.:
        raise Exception('Bad exposure.')

    # temp
    for d in gdat.indxregi:
        gdat.expo[d][where(gdat.expo[d] < 1e-50)] = 1e-50
    
    # exclude voxels with vanishing exposure
    if gdat.correxpo:
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    gdat.indxpixlrofi = intersect1d(gdat.indxpixlrofi, where(gdat.expo[d][i, :, m] > 0.)[0])
    
    gdat.indxcuberofi = meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
    
    gdat.numbpixl = gdat.indxpixlrofi.size
   
    gdat.minmexpo = empty(gdat.numbregi)
    gdat.maxmexpo = empty(gdat.numbregi)
    for d in gdat.indxregi:
        gdat.minmexpo[d] = amin(gdat.expo[d][where(gdat.expo[d] > 1e-100)])
        gdat.maxmexpo[d] = amax(gdat.expo[d])
    gdat.minmexpo = amin(gdat.minmexpo)
    gdat.maxmexpo = amax(gdat.maxmexpo)

    if gdat.minmexpo > 0:
        gdat.indxpixlroficnvt = arange(gdat.numbpixlfull)
    else:
        if gdat.verbtype > 0:
            print 'Calculating lookup table for zero-exposure pixels...'
        cntr = 0
        gdat.indxpixlroficnvt = full(gdat.numbpixlfull, -1)
        for j in gdat.indxpixlfull:
            if j in gdat.indxpixlrofi:
                gdat.indxpixlroficnvt[j] = cntr
                cntr += 1
    
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.indxcube = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
   
    gdat.indxpixleval = [[[] for d in gdat.indxregi] for y in gdat.indxgrid]
    gdat.indxcubeeval = [[[] for d in gdat.indxregi] for y in gdat.indxgrid]
    for y in gdat.indxgrid:
        for d in gdat.indxregi:
            gdat.indxpixleval[y][d] = gdat.indxpixl
            gdat.indxcubeeval[y][d] = gdat.indxcube
        
    if gdat.datatype == 'inpt':
        for d in gdat.indxregi:
            gdat.sbrtdata[d] = gdat.sbrtdata[d][gdat.indxcuberofi]

    ## exposure
    if gdat.correxpo:
        gdat.expofull = copy(gdat.expo)
        for d in gdat.indxregi:
            gdat.expo[d] = gdat.expo[d][gdat.indxcuberofi]
    
    ## backgrounds
    for strgmodl in gdat.liststrgmodl:
        indxback = getattr(gdat, strgmodl + 'indxback')
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        
        if gdat.pixltype == 'heal':
            sbrtbackhealfull = copy(sbrtbacknorm[meshgrid(indxback, gdat.indxregi, gdat.indxener, gdat.indxpixlfull, gdat.indxevtt, indexing='ij')])
            setattr(gdat, strgmodl + 'sbrtbackhealfull', sbrtbackhealfull)
        sbrtbacknormincl = [[[] for d in gdat.indxregi] for c in indxback]
        for c in indxback:
            for d in gdat.indxregi:
                sbrtbacknormincl[c][d] = sbrtbacknorm[c][d][gdat.indxcuberofi]
        setattr(gdat, strgmodl + 'sbrtbacknorm', sbrtbacknormincl)
                
    gdat.expototl = [[] for d in gdat.indxregi]
    gdat.expototlmean = [[] for d in gdat.indxregi]
    for d in gdat.indxregi:
        gdat.expototl[d] = sum(gdat.expo[d], axis=2)
        gdat.expototlmean[d] = mean(gdat.expototl[d], axis=1)

    if 'locl' in gdat.commelemspatevaltype or 'loclhash' in gdat.commelemspatevaltype:
        if gdat.exprtype == 'sdyn':
            gdat.maxmangl = 1.
        if gdat.exprtype == 'ferm':
            gdat.maxmangl = 20. / gdat.anglfact
        if gdat.exprtype == 'chan':
            gdat.maxmangl = 8. / gdat.anglfact
        if gdat.exprtype == 'hubb':
            gdat.maxmangl = 1. / gdat.anglfact
        retr_axis(gdat, 'angl', 0., gdat.maxmangl, gdat.numbangl)
    
    gdat.listnamespatmean = ['full']
    gdat.numbspatmean = len(gdat.listnamespatmean)
    gdat.indxspatmean = arange(gdat.numbspatmean)
    gdat.listindxcubespatmean = [[] for b in gdat.indxspatmean]
    for b in gdat.indxspatmean:
        if b == 0:
            gdat.listindxcubespatmean[b] = gdat.indxcube
    
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
                gdat.pixlcnvt = zeros(gdat.numbpixlfull, dtype=int) - 1
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
            numb = max(sum(gdat.fittmaxmnumbelem), sum(gdat.truemaxmnumbelem)) + 2
        else:
            numb = sum(gdat.fittmaxmnumbelem) + 2
        for k in range(numb):
            gdat.listindxpixl.append([])
            for kk in range(k):
                gdat.listindxpixl[k].append(gdat.indxpixl)
        
        # spatial averaging setup
        # temp
        if 'lghtpnts' in gdat.commelemtype or 'clus' in gdat.commelemtype:
            gdat.psfnexpr = retr_psfn(gdat, gdat.psfpexpr, gdat.indxener, gdat.binsangl, gdat.psfntypeexpr, gdat.binsoaxi, gdat.exproaxitype)
        
        # temp -- check if 1000 is too much
        gdat.numbanglelem = 1000
    
    if gdat.proplenp == None:
        if gdat.fittlensmodltype != 'none':
            gdat.proplenp = True
        else:
            gdat.proplenp = False

    # turn off relevant proposal types
    gdat.indxfixpprop = []
    for k, strg in enumerate(gdat.fittnamefixp):
        if gdat.fittnumbtrap > 0 and k in gdat.fittindxfixpnumbelemtotl:
            thisbool = False
        else:
            if strg.startswith('numbelem'):
                continue

            if gdat.fittnumbtrap > 0:
                if k in gdat.fittindxfixpmeanelem:
                    strgtemp = 'meanelem'
                if k in gdat.fittindxfixpdist:
                    strgtemp = 'dist'
            if gdat.fittpsfnevaltype != 'none':
                if k in gdat.fittindxfixppsfp:
                    strgtemp = 'psfp'
            if k in gdat.fittindxfixpbacp:
                strgtemp = 'bacp'
            if k in gdat.fittindxfixplenp:
                strgtemp = 'lenp'
            thisbool = getattr(gdat, 'prop' + strgtemp)
        
        if thisbool:
            gdat.indxfixpprop.append(getattr(gdat, 'fittindxfixp' + strg))
    gdat.propfixp = gdat.propmeanelem or gdat.propdist or gdat.propbacp or gdat.proppsfp or gdat.proplenp
    if not gdat.propfixp and not gdat.propcomp:
        raise Exception('Either a fixed dimensional parameter or an element component must be perturbed.')
    gdat.indxfixpprop = array(gdat.indxfixpprop) 
    gdat.numbfixpprop = gdat.indxfixpprop.size
    gdat.indxprop = arange(gdat.numbfixpprop)
    
    if 'lghtline' in gdat.fittelemtype:
        gdat.indxstdpelin = gdat.numbfixpprop
        gdat.indxstdpflux = gdat.numbfixpprop + 1
    else:
        gdat.indxstdplgal = gdat.numbfixpprop
        gdat.indxstdpbgal = gdat.numbfixpprop + 1
        if True in gdat.fittboolelemspec:
            gdat.indxstdpflux = gdat.numbfixpprop + 2
            if gdat.numbener > 1:
                gdat.indxstdpsind0001 = gdat.numbfixpprop + 3
                gdat.indxstdpsind0002 = gdat.numbfixpprop + 4
                gdat.indxstdpsind = gdat.numbfixpprop + 3
                gdat.indxstdpcurv = gdat.numbfixpprop + 4
                gdat.indxstdpexpc = gdat.numbfixpprop + 4
        if 'lens' in gdat.fittelemtype:
            gdat.indxstdpdefs = gdat.numbfixpprop + 2
            gdat.indxstdpasca = gdat.numbfixpprop + 3
            gdat.indxstdpacut = gdat.numbfixpprop + 4
        if 'clus' in gdat.fittelemtype:
            gdat.indxstdpnobj = gdat.numbfixpprop + 2
    # temp
    #gdat.indxstdpampl = getattr(gdat, 'indxstdp' + gdat.namefeatampl)
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
    indxelemfull = [[range(gdat.fittmaxmnumbelem[l][d]) for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl]
    indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, 'fitt')
    gdat.indxstdppara = zeros(gdat.fittnumbpara, dtype=int) - 1
    cntr = 0
    gdat.indxparaprop = zeros(gdat.numbfixpprop, dtype=int)
    for k in gdat.fittindxpara:
        if k in gdat.fittindxfixp:
            if k in gdat.indxfixpprop:
                gdat.indxstdppara[k] = cntr
                gdat.indxparaprop[cntr] = k
                cntr += 1
        else:
            for l in gdat.fittindxpopl:
                for d in gdat.fittindxregipopl[l]:
                    for strgcomp in gdat.fittliststrgcomp[l]:
                        if gdat.propcomp and k in indxsampcomp[strgcomp][l][d]:
                            gdat.indxstdppara[k] = getattr(gdat, 'indxstdp' + strgcomp)
    
    # for the fitting model, define proposal type indices
    dicttemptemp = deepcopy(gdat.__dict__) 
    for name, valu in dicttemptemp.iteritems():
        if name.startswith('fittindxfixp') and name != 'fittindxfixp' and not name.startswith('fittindxfixpnumbelem'):
            indxstdp = gdat.indxstdppara[valu]
            setattr(gdat, 'indxstdp' + name[12:], indxstdp)
    
    # proposal scale
    if 'lens' in gdat.fittelemtype:
        gdat.stdvstdp = 1e-4 + zeros(gdat.numbstdp)

        if gdat.fittnumbtrap > 0:
            gdat.stdvstdp[gdat.indxstdpmeanelempop0] = 1e-1
            gdat.stdvstdp[gdat.indxstdpdefsdistsloppop1] = 1e-1

        gdat.stdvstdp[gdat.indxstdpsigcene0evt0] = 3e-2
        gdat.stdvstdp[gdat.indxstdpbacpback0000reg0ene0] = 1e-3
        #gdat.stdvstdp[gdat.indxstdpbacpback0000reg0ene0] = 1e-1
        
        if gdat.fittlensmodltype != 'none':
            gdat.stdvstdp[gdat.indxstdplgalsour] = 1e-3
            gdat.stdvstdp[gdat.indxstdpbgalsour] = 1e-3
            gdat.stdvstdp[gdat.indxstdpfluxsour] = 1e-2
            if gdat.numbener > 1:
                gdat.stdvstdp[gdat.indxstdpsindsour] = 1e-3
            gdat.stdvstdp[gdat.indxstdpsizesour] = 1e-1
            gdat.stdvstdp[gdat.indxstdpellpsour] = 1e-1
            gdat.stdvstdp[gdat.indxstdpanglsour] = 1e-1
        if gdat.fitthostemistype != 'none':
            gdat.stdvstdp[gdat.indxstdplgalhost] = 3e-4
            gdat.stdvstdp[gdat.indxstdpbgalhost] = 3e-4
            gdat.stdvstdp[gdat.indxstdpfluxhost] = 1e-3
            if gdat.numbener > 1:
                gdat.stdvstdp[gdat.indxstdpsindhost] = 1e-3
            gdat.stdvstdp[gdat.indxstdpsizehost] = 3e-3
        if gdat.fittlensmodltype != 'none':
            gdat.stdvstdp[gdat.indxstdpbeinhost] = 1e-3
        if gdat.fitthostemistype != 'none':
            gdat.stdvstdp[gdat.indxstdpellphost] = 1e-2
            gdat.stdvstdp[gdat.indxstdpanglhost] = 1e-2
            gdat.stdvstdp[gdat.indxstdpserihost] = 1e-2
        if gdat.fittlensmodltype != 'none':
            gdat.stdvstdp[gdat.indxstdpsherextr] = 1e-1
            gdat.stdvstdp[gdat.indxstdpsangextr] = 3e-2
        
        gdat.stdvstdp[gdat.indxstdpcomp] = 5e-2
    else:
        if gdat.exprtype == 'ferm':
            gdat.stdvstdp = 1e-2 + zeros(gdat.numbstdp)
            for i in gdat.indxener:
                for c in gdat.fittindxback:
                    if c == 0:
                        gdat.stdvstdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback%04dreg%dene%d' % (c, 0, i))]] = 3e-3
                    if c == 1:
                        gdat.stdvstdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback%04dreg%dene%d' % (c, 0, i))]] = 1e-3
                    if c == 2:
                        gdat.stdvstdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpback%04dreg%dene%d' % (c, 0, i))]] = 1e-1
            if gdat.fittnumbtrap > 1:
                gdat.stdvstdp[gdat.indxstdppara[gdat.fittindxfixpmeanelem]] = 1e-2
        if gdat.exprtype == 'chan':
            gdat.stdvstdp = 1e-2 + zeros(gdat.numbstdp)
            if gdat.fittnumbtrap > 1:
                gdat.stdvstdp[gdat.indxstdpcomp] = 1e-2
                gdat.stdvstdp[gdat.indxstdpflux] = 1e-2
        if gdat.exprtype == 'sdyn':
            gdat.stdvstdp = 1e-2 + zeros(gdat.numbstdp)
    
    if gdat.sqzeprop:
        gdat.stdvstdp *= 1e-100
    if gdat.explprop:
        gdat.stdvstdp[:] = 1.

    ## input data
    if 'lghtpnts' in gdat.commelemtype or 'clus' in gdat.commelemtype:
        gdat.limsangl = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        gdat.limspsfn = [[[] for m in gdat.indxevtt] for i in gdat.indxener]
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                if gdat.trueoaxitype:
                    psfn = gdat.psfnexpr[i, :, m, 0]
                else:
                    psfn = gdat.psfnexpr[i, :, m]
                maxmpsfn = amax(psfn)
                gdat.limsangl[i][m] = [0., gdat.binsangl[amax(where(psfn > 1e-6 * maxmpsfn)[0])] * gdat.anglfact]
                gdat.limspsfn[i][m] = [maxmpsfn * 1e-6, maxmpsfn]
       
    # pixels whose posterior predicted emission will be saved
    #gdat.numbpixlsave = min(10000, gdat.numbpixl)
    #gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)
    #gdat.indxtesssave = meshgrid(gdat.indxener, gdat.indxpixlsave, gdat.indxevtt, indexing='ij')
   
    if gdat.verbtype > 1 and boolinitsetp:
        print 'fixp'
        # temp
        for strgmodl in gdat.liststrgmodl:
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
            
    if 'lghtpnts' in gdat.commelemtype or 'clus' in gdat.commelemtype:
        gdat.exprfwhm = 2. * retr_psfnwdth(gdat, gdat.psfnexpr, 0.5)
        
        # temp -- what is this?
        #gdat.stdvspatprio = amax(gdat.exprfwhm)
    #if 'lens' in gdat.commelemtype:
        # temp -- what is this?
        #gdat.stdvspatprio = amax(gdat.psfpexpr)
    
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
            # temp
            gdat.probtran = 0.4
        else:
            gdat.probtran = 0.
    gdat.probbrde = gdat.probtran - gdat.probspmr
       
    cntr = tdpy.util.cntr()
    gdat.indxproptypecomp = []
    for k in gdat.indxstdp:    
        gdat.indxproptypewith = cntr.incr()
        if k in arange(gdat.indxfixpprop.size):
            gdat.lablproptype = append(gdat.lablproptype, gdat.fittlablfixp[gdat.indxfixpprop[k]])
        else:
            for l in gdat.fittindxpopl:
                for strgcomp in gdat.fittliststrgcomp[l]:
                    if k == getattr(gdat, 'indxstdp' + strgcomp):
                        gdat.lablproptype = append(gdat.lablproptype, getattr(gdat, 'labl' + strgcomp))
            
        indx = where(gdat.indxstdppara == k)[0]
        if indx.size > 0:
            gdat.legdproptype = append(gdat.legdproptype, gdat.fittnamepara[indx[0]])
            gdat.nameproptype = append(gdat.nameproptype, gdat.namestdp[k])
        if 'lghtline' in gdat.commelemtype:
            if k >= gdat.indxstdpelin:
                gdat.indxproptypecomp.append(gdat.indxproptypewith)
        else:
            if k >= gdat.indxstdplgal:
                gdat.indxproptypecomp.append(gdat.indxproptypewith)
   
    #### filter for model elements
    gdat.listnamefilt = ['']
    if gdat.priofactdoff != 1.:
        gdat.listnamefilt += ['pars']
    #### model elements inside the image
    if 'lghtpnts' in gdat.commelemtype or 'clus' in gdat.commelemtype:
        gdat.listnamefilt += ['bndr']
    #### model subhalos inside high normalized relevance region
    if 'lens' in gdat.commelemtype:
        gdat.listnamefilt += ['nrel']
    
    if gdat.fittnumbtrap > 0.:
        
        gdat.indxproptypebrth = cntr.incr()
        gdat.indxproptypedeth = cntr.incr()
        if gdat.probbrde > 0.:
            # birth
            gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{B}$')
            gdat.legdproptype = append(gdat.legdproptype, 'Birth')
            gdat.nameproptype = append(gdat.nameproptype, 'brth')
            
            # death
            gdat.lablproptype = append(gdat.lablproptype, r'$\mathcal{D}$')
            gdat.legdproptype = append(gdat.legdproptype, 'Death')
            gdat.nameproptype = append(gdat.nameproptype, 'deth')
            
        if gdat.probspmr > 0.:
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
                    
    if gdat.datatype == 'inpt':
        proc_cntpdata(gdat)
    
    gdat.boolhash = False
    for strgmodl in gdat.liststrgmodl:
        elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        numbtrap = getattr(gdat, strgmodl + 'numbtrap')
        if numbtrap > 0:
            for l in indxpopl:
                if elemspatevaltype[l] == 'loclhash':
                    if strgmodl == 'true' and gdat.truenumbelempopl[l] > 0:
                        gdat.boolhash = True
                    if strgmodl == 'fitt' and gdat.fittmaxmnumbelempopl[l] > 0:
                        gdat.boolhash = True

    if gdat.boolhash:
        gdat.numbprox = 3
        gdat.indxprox = arange(gdat.numbprox)
        gdat.binsprox = logspace(log10(gdat.fittminmflux), log10(gdat.fittmaxmflux), gdat.numbprox + 1)
        
        # determine the maximum angle at which the contribution of the element will be computed
        gdat.maxmangleval = empty(gdat.numbprox)
        for h in gdat.indxprox:
            if gdat.specfraceval == 0:
                gdat.maxmangleval[h] = 3. * gdat.maxmgang
            else:  
                frac = min(1e-2, gdat.specfraceval * gdat.binsprox[0] / gdat.binsprox[h+1])
                psfnwdth = retr_psfnwdth(gdat, gdat.psfnexpr, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
        
        if gdat.exprtype == 'chan':
            gdat.maxmangleval = maximum(gdat.maxmangleval, array([3., 5., 7.]) / gdat.anglfact)
       
        if ('lghtpnts' in gdat.commelemtype or 'clus' in gdat.commelemtype) and gdat.maxmangl - amax(gdat.maxmangleval) < 1.1 * sqrt(2) * (gdat.maxmgang - gdat.maxmgangdata):
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
                        if gdat.maxmangl < sqrt(2.) * gdat.maxmgang:
                            raise Exception('Angular axis used to interpolate the PSF should be longer.')
                    gdat.indxpixlprox[h].append(indxpixlproxtemp)
                cntrsave = tdpy.util.show_prog(j, gdat.numbpixl, cntrsave)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
        
        gdat.numbpixlprox = zeros(gdat.numbprox) 
        for h in range(gdat.numbprox):
            for m in range(len(gdat.indxpixlprox[h])):
                if isinstance(gdat.indxpixlprox[h][m], int):
                    gdat.numbpixlprox[h] += gdat.numbpixl
                else:
                    gdat.numbpixlprox[h] += len(gdat.indxpixlprox[h][m])
            gdat.numbpixlprox[h] /= len(gdat.indxpixlprox[h])
        
        if gdat.diagmode:
            diff = gdat.numbpixlprox - roll(gdat.numbpixlprox, 1)
            if not (diff[1:] >= 0).all():
                print 'gdat.numbpixlprox'
                print gdat.numbpixlprox
                print 'roll(gdat.numbpixlprox, 1)'
                print roll(gdat.numbpixlprox, 1)
                print 'diff'
                print diff
                raise Exception('Number of pixels in the pixel look-up table is wrong.')
    

def rplc_list(listarry, indx):
    
    listarryrplc = [[] for k in range(len(listarry))]
    for k in range(len(listarry)):
        listarryrplc[k] = listarry[k][indx]

    return listarryrplc


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
                            temp = retr_rele(gdat, gdat.truecntp['lens'][0, :, 0], gdat.binslgal[k], gdat.binsbgal[n], 1., \
                                                                                                    gdat.trueascadistmeanpop0, gdat.trueacutdistmeanpop0, gdat.indxpixl)
    
                            #temp /= amax(temp)
                            pdfnpriotemp[k, n] = temp**gdat.relnpowr
                    lpdfprio, lpdfprioobjt = retr_spatprio(gdat, pdfnpriotemp)
                    lpdfpriointp = lpdfprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
        setattr(gdat, strgmodl + 'lpdfprio', lpdfprio)
        setattr(gdat, strgmodl + 'lpdfprioobjt', lpdfprioobjt)
        setattr(gdat, strgmodl + 'lpdfpriointp', lpdfpriointp)
        
    # find the pixels at which data count maps have local maxima
    if gdat.pixltype == 'cart':
        for d in gdat.indxregi:
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    # temp
                    gdat.indxxdatmaxm, gdat.indxydatmaxm = tdpy.util.retr_indximagmaxm(gdat.cntpdatacart[d, i, :, m])

    # sanity checks
    # temp
    if (fabs(gdat.cntpdata - rint(gdat.cntpdata)) > 1e-3).any() and boolinitsetp:
        print 'Fractional counts!'

    if amin(gdat.cntpdata) < 0. and boolinitsetp:
        print 'Negative counts!'

    if gdat.jitt:
        gdat.llik = empty_like(gdat.cntpdata)
        cntpmodl = zeros_like(gdat.cntpdata)
        gdat.numbthre = 4
        gdat.indxthre = arange(gdat.numbthre)
        gdat.sizechun = (gdat.numbpixl + gdat.numbthre - 1) // gdat.numbthre
        gdat.llikchun = [gdat.llik[:, k*gdat.sizechun:(k+1)*gdat.sizechun, :] for k in gdat.indxthre]
        gdat.cntpdatachun = [gdat.cntpdata[:, k*gdat.sizechun:(k+1)*gdat.sizechun, :] for k in gdat.indxthre]
        timefunc(None, "numpy (1 thread)", retr_llik_nump, gdat, cntpmodl)
        timefunc("numba (1 thread)", retr_llik_sequ, gdat.llik, gdat.cntpdata, cntpmodl)
        timefunc("numba (%d threads)" % gdat.numbthre, retr_llik_mult, gdat, cntpmodl)


def retr_gradmaps(gdat, maps):
    
    # temp -- this does not work with vanishing exposure
    maps = maps.reshape((gdat.numbsidecart, gdat.numbsidecart))
    grad = dstack((gradient(maps, gdat.sizepixl, axis=0), gradient(maps, gdat.sizepixl, axis=1))).reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
    grad = grad.reshape((gdat.numbpixl, 2))

    return grad


def retr_spatmean(gdat, inpt, indxregieval, boolcntp=False):
    
    listspatmean = empty((gdat.numbspatmean, gdat.numbregi, gdat.numbener))
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        for dd, d in enumerate(indxregieval):
            if boolcntp:
                cntp = inpt[dd][gdat.listindxcubespatmean[b]]
            else:
                cntp = inpt[dd][gdat.listindxcubespatmean[b]] * gdat.expo[d][gdat.listindxcubespatmean[b]] * gdat.apix
                if gdat.enerdiff:
                    cntp *= gdat.deltener[:, None, None]
            
            spatmean = mean(sum(cntp, 2), axis=1) / gdat.expototlmean[d] / gdat.apix
            if gdat.enerdiff:
                spatmean /= gdat.deltener
            
            listspatmean[b, :, :] = spatmean

    return listspatmean


def retr_rele(gdat, maps, lgal, bgal, defs, asca, acut, indxpixleval, absv=True):
    
    grad = retr_gradmaps(gdat, maps)
        
    defl = retr_defl(gdat, lgal, bgal, defs, 0., 0., asca=asca, acut=acut, indxpixleval=indxpixleval)
    dotstemp = sum(grad * defl, 1)
    if absv:
        dotstemp = fabs(dotstemp)
    else:
        dotstemp = dotstemp
    
    dots = mean(dotstemp) * gdat.sizepixl
    
    return dots


def retr_llik_nump(gdat, cntpmodl):
    
    llik = gdat.cntpdata * log(cntpmodl) - cntpmodl
    
    return llik


@jit(nopython=True, nogil=True)
def retr_llik_sequ(llik, cntpdata, cntpmodl):
    
    for i in range(llik.shape[0]):
        for j in range(llik.shape[1]):
            for m in range(llik.shape[2]):
                llik[i, j, m] = cntpdata[i, j, m] * log(cntpmodl[i, j, m]) - cntpmodl[i, j, m]
    

def retr_llik_mult(gdat, cntpmodl):
    
    listchun = []
    for k in gdat.indxthre:
        listchun.append([gdat.llikchun[k], gdat.cntpdatachun[k], cntpmodl[:, k*gdat.sizechun:(k+1)*gdat.sizechun, :]])
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
    gdat.maxmcntpdata = max(0.5, amax(sum(gdat.cntpdata, 3)))
    gdat.maxmcntpmodl = gdat.maxmcntpdata
    gdat.minmcntpmodl = 1e-3 * gdat.maxmcntpdata
    gdat.minmcntpdata = max(0.5, gdat.minmcntpmodl)
    gdat.maxmcntpresi = ceil(gdat.maxmcntpdata * 0.1)
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
        minmscal = arcsinh(minmscal)
    if scal == 'logt':
        minmscal = log10(minmscal)
    maxmscal = maxm
    if scal == 'asnh':
        maxmscal = arcsinh(maxmscal)
    if scal == 'logt':
        maxmscal = log10(maxmscal)

    tickscal = linspace(minmscal, maxmscal, gdat.numbtickcbar)
    labl = empty(gdat.numbtickcbar, dtype=object)
    tick = copy(tickscal)
    for k in range(gdat.numbtickcbar):
        if scal == 'asnh':
            tick[k] = sinh(tickscal[k])
        elif scal == 'logt':
            tick[k] = 10**(tickscal[k])

        # avoid very small, but nonzero central values in the residual count color maps
        if strgcbar == 'cntpresi' and fabs(tick[k]) < 1e-5:
            tick[k] = 0.

        if strgcbar == 'cntpdata' and amax(tick) > 1e3:
            labl[k] = '%d' % tick[k]
        else:
            labl[k] = '%.3g' % tick[k]
    
    setattr(gdat, 'tick' + strgcbar, tickscal)
    setattr(gdat, 'minmscal' + strgcbar, minmscal)
    setattr(gdat, 'maxmscal' + strgcbar, maxmscal)
    setattr(gdat, 'labl' + strgcbar, labl)
    

def retr_fromgdat(gdat, gdatmodi, strgstat, strgmodl, strgvarb, mometype='medi', indxvarb=None, indxlist=None):
    
    if strgvarb == 'cntpdata' or strgvarb.startswith('histcntpdata'):
        varb = getattr(gdat, strgvarb)
    else:
        if strgmodl == 'true':
            varb = getattr(gdat, 'true' + strgvarb)
        if strgmodl == 'fitt':
            if strgstat == 'this':
                if mometype == 'errr':
                    varb = getattr(gdatmodi, strgstat + 'errr' + strgvarb)
                else:
                    varb = getattr(gdatmodi, strgstat + strgvarb)
            if strgstat == 'post':
                varb = getattr(gdat, mometype + strgvarb)
    
    if indxlist != None:
        varb = varb[indxlist]

    if indxvarb != None:
        if mometype == 'errr':
            varb = varb[[slice(None)] + indxvarb]
        else:
            varb = varb[indxvarb]

    return copy(varb)


def retr_indxsamp(gdat, strgmodl='fitt', init=False):
    
    dicttemp = {}
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    if init:
        # transdimensional element populations
        backtype = getattr(gdat, strgmodl + 'backtype')
        
        for c in range(len(backtype)):
            if isinstance(backtype[c], str):
                if backtype[c].startswith('bfunfour') or backtype[c].startswith('bfunwfou'):
                    namebfun = backtype[c][:8]
                    ordrexpa = int(backtype[c][8:])
                    numbexpa = 4 * ordrexpa**2
                    indxexpa = arange(numbexpa)
                    del backtype[c]
                    for k in indxexpa:
                        backtype.insert(c+k, namebfun + '%04d' % k)
        numbback = len(backtype)
        numbpopl = len(elemtype)
        print 'numbback'
        print numbback
        print 'backtype'
        print backtype
        
        indxpopl = arange(numbpopl)
        indxback = arange(numbback)
        
        # feature used to model the amplitude of elements
        namefeatampl = [[] for l in indxpopl]
        indxcompampl = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l] == 'lghtline':
                namefeatampl[l] = 'flux'
                indxcompampl[l] = 1
            elif elemtype[l].startswith('lght'):
                namefeatampl[l] = 'flux'
                indxcompampl[l] = 2
            if elemtype[l] == 'lens':
                namefeatampl[l] = 'defs'
                indxcompampl[l] = 2
            if elemtype[l] == 'clus':
                namefeatampl[l] = 'nobj'
                indxcompampl[l] = 2
    else:
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        elemspatevaltype = getattr(gdat, strgmodl + 'elemspatevaltype')
        namefeatampl = getattr(gdat, strgmodl + 'namefeatampl') 
        indxcompampl = getattr(gdat, strgmodl + 'indxcompampl') 
        maxmnumbelem = getattr(gdat, strgmodl + 'maxmnumbelem') 
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        indxback = getattr(gdat, strgmodl + 'indxback')
        numbregipopl = getattr(gdat, strgmodl + 'numbregipopl')
        indxregipopl = getattr(gdat, strgmodl + 'indxregipopl')
        hostemistype = getattr(gdat, strgmodl + 'hostemistype')
        listnameback = getattr(gdat, strgmodl + 'listnameback') 
        
        if strgmodl == 'true':
            numbelempopl = zeros(numbpopl)
            for l in indxpopl:
                for d in indxregipopl[l]:
                    numbelempopl[l] += getattr(gdat, 'truenumbelempop%dreg%d' % (l, d))
         
        numbelemzero = [[] for l in indxpopl]
        for l in indxpopl:
            numbelemzero[l] = zeros(numbregipopl[l], dtype=int)

        # element setup
        ## flag to calculate the kernel approximation errors
        calcerrr = [[] for l in indxpopl]
        if elemspatevaltype[l] != 'full' and gdat.numbpixlfull < 1e5:
            # temp
            #calcerrr[l] = True
            calcerrr[l] = False
        else:
            calcerrr[l] = False
        
        ## sorting feature
        namefeatsort = [[] for l in indxpopl]
        for l in indxpopl:
            # feature to be used to sort elements
            if elemtype[l].startswith('lght'):
                namefeatsort[l] = 'flux'
            if elemtype[l] == 'lens':
                namefeatsort[l] = 'defs'
            if elemtype[l] == 'clus':
                namefeatsort[l] = 'nobj'
    
        ## selection feature
        listnamefeatsele = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght'):
                listnamefeatsele[l] = ['flux']
            if elemtype[l] == 'lens':
                listnamefeatsele[l] = ['defs', 'mcut', 'rele']
            if elemtype[l] == 'clus':
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
                if elemtype[l] == 'lghtpnts':
                    lablelemextn[l] = r'\rm{pts}'
            if elemtype == 'lens':
                lablelemextn[l] = r'\rm{sub}'
            if elemtype == 'clus':
                lablelemextn[l] = r'\rm{cls}'
            if elemtype == 'lghtline':
                lablelemextn[l] = r'\rm{lin}'
    
        ## legends
        listlegdelem = [[] for l in indxpopl]
        for l in indxpopl:
            if gdat.numbgrid > 1:
                if elemtype[l] == 'lghtpnts':
                    listlegdelem[l] = 'FPS'
                if elemtype[l] == 'lghtgausbgrd':
                    listlegdelem[l] = 'BPS'
            else:
                if elemtype[l] == 'lghtpnts':
                    listlegdelem[l] = 'PS'
            if elemtype == 'lens':
                listlegdelem[l] = 'Subhalo'
            if elemtype == 'clus':
                listlegdelem[l] = 'Cluster'
            if elemtype[l] == 'lghtline':
                listlegdelem[l]= 'Line'
    
        indxpoplgrid = [[] for y in gdat.indxgrid]
        for y in gdat.indxgrid: 
            for indx, elemtypetemp in enumerate(elemtype):
                # foreground grid (image plane) -- the one where the data is measured
                if y == 0:
                    if elemtypetemp.startswith('lght') and not elemtypetemp.endswith('bgrd') or elemtypetemp == 'clus':
                        indxpoplgrid[y].append(indx)
                # foreground mass grid
                if y == 1:
                    if elemtypetemp.startswith('lens'):
                        indxpoplgrid[y].append(indx)
                # background grid (source plane)
                if y == 2:
                    if elemtypetemp.endswith('bgrd')
                        indxpoplgrid[y].append(indx)
        
        elemregitype = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('hypr'):
                elemregitype[l] = False
            else:
                elemregitype[l] = True
        
        indxgridpopl = [[] for l in indxpopl]
        for l in indxpopl:
            for y in gdat.indxgrid:
                if l in indxpoplgrid[y]:
                    indxgridpopl[l] = y
    
        if 'lghtpnts' in elemtype:
            calcelemsbrt = True
        else:
            calcelemsbrt = False

        if 'lght' in elemtype:
            calcelemsbrtbgrd = True
        else:
            calcelemsbrtbgrd = False

        if 'lens' in elemtype:
            calcelemdefl = True
        else:
            calcelemdefl = False

        ## element Boolean flags
        boolelemspec = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght'):
                boolelemspec[l] = True
            else:
                boolelemspec[l] = False
        boolelemspectotl = True in boolelemspec

        boolelemfore = [[] for l in indxpopl]
        for l in indxpopl:
            if boolelemspec[l] or elemtype[l] == 'clus':
                boolelemfore[l] = True
            else:
                boolelemfore[l] = False
    
        boolelempsfn = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l] == 'lghtpnts' or elemtype[l] == 'clus':
                boolelempsfn[l] = True
            else:
                boolelempsfn[l] = False
        
        spectype = [[] for l in indxpopl]
        for l in indxpopl:
            if boolelemspec[l]:
                spectype[l] = 'powr'
            else:
                spectype[l] = 'none'
        setp_varbvalu(gdat, 'spectype', spectype, strgmodl=strgmodl)
    
        if boolelemspectotl:
            # flux
            if gdat.exprtype == 'ferm':
                minmflux = 3e-9
                maxmflux = 1e-6
            if gdat.exprtype == 'chan':
                if gdat.anlytype == 'spec':
                    minmflux = 1e-7
                else:
                    minmflux = 3e-9
                maxmflux = 1e-6
            if gdat.exprtype == 'sdyn':
                minmflux = 0.1
                maxmflux = 100.
            if gdat.exprtype == 'hubb':
                minmflux = 0.1
                maxmflux = 100.
            setp_varblimt(gdat, 'minmflux', [minmflux, maxmflux], strgmodl=strgmodl)
    
            ### spectral parameters
            if gdat.exprtype == 'ferm':
                sind = [1., 3.]
                minmsind = 1.
                maxmsind = 3.
                minmsind0001 = -5.
                maxmsind0001 = 10.
                minmsind0002 = -5.
                maxmsind0002 = 10.
                setp_varblimt(gdat, 'sind0001', [minmsind0001, maxmsind0001], strgmodl=strgmodl)
                setp_varblimt(gdat, 'sind0002', [minmsind0002, maxmsind0002], strgmodl=strgmodl)
            if gdat.exprtype == 'chan':
                minmsind = 0.4
                maxmsind = 2.4
                sind = [0.4, 2.4]
            if gdat.exprtype == 'hubb':
                minmsind = 0.5
                maxmsind = 2.5
                sind = [0.4, 2.4]
            setp_varblimt(gdat, 'sind', [minmsind, maxmsind], strgmodl=strgmodl)
            setp_varblimt(gdat, 'sinddistmean', sind, popl='full', strgmodl=strgmodl)
            #### standard deviations should not be too small
            setp_varblimt(gdat, 'sinddiststdv', [0.3, 2.], popl='full', strgmodl=strgmodl)
            setp_varblimt(gdat, 'curvdistmean', [-1., 1.], popl='full', strgmodl=strgmodl)
            setp_varblimt(gdat, 'curvdiststdv', [0.1, 1.], popl='full', strgmodl=strgmodl)
            setp_varblimt(gdat, 'expcdistmean', [1., 8.], popl='full', strgmodl=strgmodl)
            setp_varblimt(gdat, 'expcdiststdv', [0.01 * gdat.maxmener, gdat.maxmener], popl='full', strgmodl=strgmodl)
    
        # total maximum number of elements
        maxmnumbelempopl = zeros(numbpopl)
        for l in indxpopl:
            maxmnumbelempopl += sum(maxmnumbelem[l])
        maxmnumbelemtotl = sum(maxmnumbelempopl) 

        # construct background surface brightness templates from the user input
        lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
        backtype = getattr(gdat, strgmodl + 'backtype')
        numbback = getattr(gdat, strgmodl + 'numbback')
        sbrtbacknorm = empty((numbback, gdat.numbregi, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
        unifback = ones(numbback, dtype=bool)
        
        for c in indxback:
            if backtype[c] == 'data':
                sbrtbacknormtemp = copy(gdat.sbrtdata)
                sbrtbacknormtemp[where(sbrtbacknormtemp == 0.)] = 1e-100
            elif isinstance(backtype[c], float):
                if gdat.pixltype == 'heal':
                    sbrtbacknormtemp = zeros((gdat.numbregi, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + backtype[c]
                if gdat.pixltype == 'cart':
                    sbrtbacknormtemp = zeros((gdat.numbregi, gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + backtype[c]
            elif isinstance(backtype[c], list):
                sbrtbacknormtemp = retr_spec(gdat, array([backtype[c][0]]), sind=array([backtype[c][1]]))[:, 0][None, :, None, None] * \
                                                                    ones((gdat.numbregi, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            elif isinstance(backtype[c], ndarray) and backtype[c].ndim == 1:
                if gdat.pixltype == 'heal':
                    sbrtbacknormtemp = zeros((gdat.numbregi, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + backtype[c][:, None, None]
                if gdat.pixltype == 'cart':
                    sbrtbacknormtemp = zeros((gdat.numbregi, gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + backtype[c][:, None, None]
            elif backtype[c].startswith('bfunfour') or backtype[c].startswith('bfunwfou'):
                namebfun = getattr(gdat, strgmodl + 'namebfun')
                ordrexpa = getattr(gdat, strgmodl + 'ordrexpa')
                numbexpa = getattr(gdat, strgmodl + 'numbexpa')
                indxexpatemp = int(backtype[c][8:]) 
                indxterm = indxexpatemp // ordrexpa**2
                indxexpaxdat = (indxexpatemp % ordrexpa**2) // ordrexpa
                indxexpaydat = (indxexpatemp % ordrexpa**2) % ordrexpa
                if namebfun == 'bfunfour':
                    ampl = 1.
                    func = gdat.meanbgalcart 
                if namebfun == 'bfunwfou':
                    ampl = sqrt(gdat.fdfmderiintp)
                    func = gdat.fdfmintp
                argslgal = 2. * pi * indxexpaxdat * gdat.meanlgalcart / gdat.maxmgangdata
                if indxterm == 0:
                    termfrst = sin(argslgal)
                    termseco = ampl * sin(2. * pi * indxexpaydat * func / gdat.maxmgangdata)
                if indxterm == 1:
                    termfrst = sin(argslgal)
                    termseco = ampl * cos(2. * pi * indxexpaydat * func / gdat.maxmgangdata)
                if indxterm == 2:
                    termfrst = cos(argslgal)
                    termseco = ampl * sin(2. * pi * indxexpaydat * func / gdat.maxmgangdata)
                if indxterm == 3:
                    termfrst = cos(argslgal)
                    termseco = ampl * cos(2. * pi * indxexpaydat * func / gdat.maxmgangdata)
                sbrtbacknormtemp = termfrst[None, :] * termseco[:, None]
                
                if False:
                    print 'backtype'
                    print backtype
                    print 'gdat.meanlgalcart'
                    summgene(gdat.meanlgalcart)
                    print 'gdat.meanbgalcart'
                    summgene(gdat.meanbgalcart)
                    print 'termfrst'
                    summgene(termfrst)
                    print 'termseco'
                    summgene(termseco)
                    print 'termfrst[None, :] * termseco[: None]'
                    summgene(termfrst[None, :] * termseco[: None])
                    print 'sbrtbacknormtemp'
                    summgene(sbrtbacknormtemp)
                    print 'numbback'
                    print numbback
                    print 
                sbrtbacknormtemptemp = empty((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
                for i in gdat.indxenerfull:
                    for m in gdat.indxevttfull:
                        sbrtbacknormtemptemp[i, :, m] = sbrtbacknormtemp.flatten()
                sbrtbacknormtemp = sbrtbacknormtemptemp
            else:
                path = gdat.pathinpt + backtype[c]
                sbrtbacknormtemp = pf.getdata(path)
                
                # temp
                if (gdat.pixltype == 'heal' or gdat.pixltype == 'cart' and gdat.forccart) and sbrtbacknormtemp.ndim == 3 or gdat.pixltype == 'cart' and sbrtbacknormtemp.ndim == 4:
                    print 'Input background template incompatible with PCAT %s. Converting...' % gdat.strgvers
                    sbrtbacknormtemp = sbrtbacknormtemp[None, ...]

                # temp
                if gdat.exprtype == 'ferm' and sbrtbacknormtemp.shape[3] == 4 and (gdat.anlytype.startswith('rec7') or gdat.anlytype.startswith('manu')):
                    print 'Input background template incompatible with new Fermi-LAT event type binning. Converting...'
                    sbrtbacknormtemp = sbrtbacknormtemp[:, :, :, 2:4]

                if gdat.pixltype == 'cart':
                    if gdat.forccart:
                        sbrtbacknormtemptemp = empty((gdat.numbregi, gdat.numbenerfull, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevttfull))
                        for d in gdat.indxregi:
                            for i in gdat.indxenerfull:
                                for m in gdat.indxevttfull:
                                    sbrtbacknormtemptemp[d, i, :, :, m] = tdpy.util.retr_cart(sbrtbacknormtemp[d, i, :, m], \
                                                                                        numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                                                        minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                                        minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).T
                        sbrtbacknormtemp = sbrtbacknormtemptemp
                    
                    if sbrtbacknormtemp.shape[2] != gdat.numbsidecart:
                        print 'gdat.numbsidecart'
                        print gdat.numbsidecart
                        print 'sbrtbacknormtemp.shape[2]'
                        print sbrtbacknormtemp.shape[2]
                        raise Exception('Provided background template must have the chosen image dimensions.')
            
                    sbrtbacknormtemp = sbrtbacknormtemp.reshape((sbrtbacknormtemp.shape[0], sbrtbacknormtemp.shape[1], -1, sbrtbacknormtemp.shape[-1]))
           
            sbrtbacknorm[c, ...] = sbrtbacknormtemp
           
            
            # determine spatially uniform background templates
            for d in gdat.indxregi:
                for i in gdat.indxenerfull:
                    for m in gdat.indxevttfull:
                        if std(sbrtbacknorm[c, d, i, :, m]) > 1e-6:
                            unifback[c] = False

        for c in indxback:
            if amin(sbrtbacknorm[c, ...]) < 0. and isinstance(backtype[c], str) and not backtype[c].startswith('bfun'):
                booltemp = False
                raise Exception('Background templates must be positive-definite everywhere.')
        
        if not isfinite(sbrtbacknorm).all():
            raise Exception('')

        boolzero = True
        for c in indxback:
            if amin(sbrtbacknorm[c, ...]) > 0. or backtype[c] == 'data':
                boolzero = False
        
        if boolzero and isinstance(backtype[0], str) and not backtype[0].startswith('bfun'):
            raise Exception('At least one background template must be positive everywhere.')
       
        if 'data' in backtype and len(backtype) != 1:
            raise Exception('data - PS residual can be the only background.')
        if backtype[0] != 'data' and gdat.penalpridiff:
            raise Exception('Diffuse background power spectrum penalization is unncessary if the background is not the data - PS residual.')
        
        if maxmnumbelemtotl > 0 and ('lghtpnts' in elemtype or 'clus' in elemtype):
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
            oaxitype = getattr(gdat, strgmodl + 'oaxitype')
        
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
            
            if gdat.psfninfoprio:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        meansigc = gdat.psfpexpr[i * numbpsfptotl + m * numbpsfptotl * gdat.numbener]
                        stdvsigc = meansigc * 0.1
                        setp_varblimt(gdat, 'sigcene%devt%d' % (i, m), [meansigc, stdvsigc], typelimt='meanstdv', strgmodl=strgmodl)
                        if psfntype == 'doubking' or psfntype == 'singking':
                            meangamc = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 1]
                            stdvgamc = meangamc * 0.1
                            setp_varblimt(gdat, 'gamcene%devt%d' % (i, m), [meangamc, stdvgamc], typelimt='meanstdv', strgmodl=strgmodl)
                            if psfntype == 'doubking':
                                meansigt = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 2]
                                stdvsigt = meansigt * 0.1
                                setp_varblimt(gdat, 'sigtene%devt%d' % (i, m), [meansigt, stdvsigt], typelimt='meanstdv', strgmodl=strgmodl)
                                meangamt = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 3]
                                stdvgamt = meangamt * 0.1
                                setp_varblimt(gdat, 'gamtene%devt%d' % (i, m), [meangamt, stdvgamt], typelimt='meanstdv', strgmodl=strgmodl)
                                meanpsff = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 4]
                                stdvpsff = meanpsff * 0.1
                                setp_varblimt(gdat, 'psffene%devt%d' % (i, m), [meanpsff, stdvpsff], typelimt='meanstdv', strgmodl=strgmodl)
                        elif oaxitype:
                            meanonor = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 1]
                            stdvonor = meanonor * 0.1
                            setp_varblimt(gdat, 'onorene%devt%d' % (i, m), [meanonor, stdvonor], typelimt='meanstdv', strgmodl=strgmodl)
                            meanoind = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 2]
                            stdvoind = meanoind * 0.1
                            setp_varblimt(gdat, 'oindene%devt%d' % (i, m), [meanoind, stdvoind], typelimt='meanstdv', strgmodl=strgmodl)
            else:
                if gdat.exprtype == 'sdyn':
                    minmsigm = 0.1 / gdat.anglfact
                    maxmsigm = 0.2 / gdat.anglfact
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
                minmonor = 0.01
                maxmonor = 1.
                minmoind = 1.5
                maxmoind = 2.5
                setp_varblimt(gdat, 'sigc', [minmsigm, maxmsigm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'sigt', [minmsigm, maxmsigm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'gamc', [minmgamm, maxmgamm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'gamt', [minmgamm, maxmgamm], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'onor', [minmonor, maxmonor], ener='full', evtt='full', strgmodl=strgmodl)
                setp_varblimt(gdat, 'oind', [minmoind, maxmoind], ener='full', evtt='full', strgmodl=strgmodl)
            setp_varblimt(gdat, 'psff', [0., 1.], ener='full', evtt='full', strgmodl=strgmodl)
 
        specback = getattr(gdat, strgmodl + 'specback')
        
        spectype = getattr(gdat, strgmodl + 'spectype')

        # background
        ## number of background parameters
        numbbacp = 0
        for c in indxback:
            if specback[c]:
                numbbacp += gdat.numbregi
            else:
                numbbacp += gdat.numbener * gdat.numbregi
   
        ## background parameter indices
        indxbackbacp = zeros(numbbacp, dtype=int)
        indxenerbacp = zeros(numbbacp, dtype=int)
        indxregibacp = zeros(numbbacp, dtype=int)
        cntr = 0
        for c in indxback:
            if specback[c]:
                for d in gdat.indxregi: 
                    indxbackbacp[cntr] = c
                    indxregibacp[cntr] = d
                    cntr += 1
            else:
                for d in gdat.indxregi: 
                    for i in gdat.indxener:
                        indxenerbacp[cntr] = i
                        indxbackbacp[cntr] = c
                        indxregibacp[cntr] = d
                        cntr += 1
        
        indxbacpback = [[] for c in indxback]
        for c in indxback:
            if specback[c]:
                indxbacpback[c] = where(indxbackbacp == c)[0]
            else:
                indxbacpback[c] = where(indxbackbacp == c)[0].reshape((gdat.numbregi, gdat.numbener))
                
        # features which correlate with significance
        liststrgfeatsign = [[] for l in indxpopl]
        for l in indxpopl:
            for namefeat in listnamefeatsele[l]:
                liststrgfeatsign[l] += [namefeat]
            
        liststrgcomp = [[] for l in indxpopl]
        listscalcomp = [[] for l in indxpopl]
        listscalcomp = [[] for l in indxpopl]
        for l in indxpopl:
            
            if elemtype[l] == 'lghtline':
                liststrgcomp[l] = ['elin']
                listscalcomp[l] = ['self']
            else:
                liststrgcomp[l] = ['lgal', 'bgal']
                listscalcomp[l] = ['self', 'self']
            
            if elemtype[l].startswith('lght'):
                liststrgcomp[l] += ['flux']
            if elemtype[l] == 'lens':
                liststrgcomp[l] += ['defs']
            if elemtype[l] == 'clus':
                liststrgcomp[l] += ['nobj']
            listscalcomp[l] += ['powrslop']
            
            if gdat.numbener > 1 and elemtype[l].startswith('lght'):
                if spectype[l] == 'colr':
                    for i in gdat.indxener:
                        if i == 0:
                            continue
                        liststrgcomp[l] += ['sind%04d' % i]
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
        
        # variables for which whose marginal distribution and pair-correlations will be plotted
        liststrgfeatodim = [[] for l in indxpopl]
        for l in indxpopl:
            liststrgfeatodim[l] = deepcopy(liststrgcomp[l])
            liststrgfeatodim[l] += ['deltllik']
            if gdat.numbpixlfull > 1:
                liststrgfeatodim[l] += ['gang', 'aang']
            if elemtype[l] == 'lght':
                liststrgfeatodim[l] += ['cnts']
                if gdat.exprtype == 'ferm':
                    liststrgfeatodim[l] + ['sbrt0018']
            if elemtype[l] == 'lens':
                liststrgfeatodim[l] += ['mcut', 'diss', 'rele', 'reln', 'reld', 'relc']

        # add reference element features that are not available in the PCAT element model
        if strgmodl == 'fitt':
            gdat.refrliststrgfeatonly = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
            for q in gdat.indxrefr: 
                for name in gdat.refrliststrgfeat[q]:
                    for l in indxpopl:
                        if not name in liststrgfeatodim[l]:
                            liststrgfeatodim[l].append(name)
                            gdat.refrliststrgfeatonly[q][l].append(name)

        # defaults
        liststrgpdfnmodu = [[] for l in indxpopl]
        liststrgfeatmodu = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l] == 'lght': 
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
        
        # temp
        for l in indxpopl:
            if elemtype[l] == 'lghtline':
                liststrgfeat[l] += ['spec']
            if elemtype[l] == 'lght':
                liststrgfeat[l] += ['spec', 'specplot']
            if elemtype[l] == 'lens':
                liststrgfeat[l] += ['deflprof']
        
        liststrgfeateval = [[] for l in indxpopl]
        for l in indxpopl:
            if elemtype[l].startswith('lght'):
                liststrgfeateval[l] = ['lgal', 'bgal', 'spec']
            if elemtype[l] == 'clus':
                liststrgfeateval[l] = ['lgal', 'bgal', 'nobj']
            if elemtype[l] == 'lens':
                liststrgfeateval[l] = ['lgal', 'bgal', 'defs', 'asca', 'acut']
            if elemtype[l] == 'lghtline':
                liststrgfeateval[l] = ['elin', 'spec']

        # variables for which pair-correlations will be plotted
        liststrgfeatcorr = [[] for l in indxpopl]
        if gdat.plotelemcorr:
            for l in indxpopl:
                for strgfeat in liststrgfeatodim[l]:
                    liststrgfeatcorr[l].append(strgfeat)
        
        # number of element parameters
        numbcomp = zeros(numbpopl, dtype=int)
        for l in indxpopl:
            numbcomp[l] = len(liststrgcomp[l])
        maxmnumbcomp = amax(numbcomp)
        
        indxcomp = []
        for l in indxpopl:
            indxcomp.append(arange(numbcomp[l]))

        # number of transdimensional parameters
        numbtrapregipopl = [[] for l in indxpopl]
        numbtrapregipoplcumr = [[] for l in indxpopl]
        for l in indxpopl:
            numbtrapregipopl[l] = maxmnumbelem[l] * numbcomp[l]
            numbtrapregipoplcumr[l] = cumsum(numbtrapregipopl[l])
            for ll in range(l):
                numbtrapregipoplcumr[l] += sum(numbtrapregipopl[ll])
        numbtrapregipoplcuml = [[] for l in indxpopl]
        for l in indxpopl:
            numbtrapregipoplcuml[l] = copy(numbtrapregipoplcumr[l])
        numbtrappoplcuml = zeros(numbpopl, dtype=int)
        numbtrappoplcumr = zeros(numbpopl, dtype=int)
        for l in indxpopl:
            numbtrappoplcuml[l] = numbtrapregipoplcuml[l][0]
            numbtrappoplcumr[l] = numbtrapregipoplcumr[l][-1]
        for l in indxpopl[::-1]:
            for d in indxregipopl[l][::-1]:
                if d == 0:
                    numbtrapregipoplcuml[l][d] = numbtrapregipoplcuml[l-1][-1]
                else:
                    numbtrapregipoplcuml[l][d] = numbtrapregipoplcuml[l][d-1]
        numbtrapregipoplcuml[0][0] = 0
        
        numbtrap = 0
        for l in indxpopl:
            numbtrap += sum(numbtrapregipopl[l])
        
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
        
        numbdeflsubhplot = 2
        numbdeflsingplot = numbdeflsubhplot
        if numbtrap > 0:
            numbdeflsingplot += 3

        listnamediff = []
        for c in indxback:
            listnamediff += ['back%04d' % c]
        if hostemistype != 'none':
            listnamediff += ['host']
        if lensmodltype != 'none':
            listnamediff += ['lens']
        
        listnameecom = deepcopy(listnamediff)
        if 'lght' in elemtype:
            listindxpopllght = [indx for indx, elemtypetemp in enumerate(elemtype) if elemtypetemp == 'lght']
            for indxpopllght in listindxpopllght:
                if strgmodl == 'true' and numbelempopl[l] > 0 or strgmodl == 'fitt' and maxmnumbelempopl[l] > 0:
                    if not 'pnts' in listnameecom:
                        listnameecom += ['pnts']
                    if not 'pntssubt' in listnameecom:
                        listnameecom += ['pntssubt']
        listnameecomtotl = listnameecom + ['modl']

        numbdiff = len(listnamediff)
        convdiff = zeros(numbdiff, dtype=bool)
        for k, namediff in enumerate(listnamediff):
            if namediff.startswith('back'):
                indx = int(namediff[-4:])
                convdiff[k] = not unifback[indx]
            else:
                convdiff[k] = True

        cntr = tdpy.util.cntr()
        
        liststrgfeatdefa = deepcopy(liststrgfeattotl)
        if 'lght' in elemtype:
            for strgfeat in ['sind', 'curv', 'expc', 'sind0001', 'sind0002']:
                if not strgfeat in liststrgfeatdefa:
                    liststrgfeatdefa.append(strgfeat)

        if numbtrap > 0:

            # number of elements
            for l in indxpopl:
                for d in indxregipopl[l]:
                    dicttemp['indxfixpnumbelempop%dreg%d' % (l, d)] = cntr.incr()
            
            # hyperparameters
            ## mean number of elements
            for l in indxpopl:
                dicttemp['indxfixpmeanelempop%d' % l] = cntr.incr()
            
            for strgfeat in liststrgfeatpriototl:
                try:
                    disttype = getattr(gdat, strgmodl + strgfeat + 'disttype')
                except:
                    if strgfeat == 'lgal' or strgfeat == 'bgal' or strgfeat == 'elin':
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
        
            dicttemp['indxfixpnumbelem'] = [[] for l in indxpopl]
            for l in indxpopl:
                dicttemp['indxfixpnumbelem'][l] = empty(numbregipopl[l], dtype=int)
            dicttemp['indxfixpmeanelem'] = []
            for strg, valu in dicttemp.iteritems():
                if strg[8:].startswith('numbelemp'):
                    indxpopltemp = int(strg[8:][11])
                    indxregitemp = int(strg[8:][15])
                    dicttemp['indxfixpnumbelem'][indxpopltemp][indxregitemp] = valu
                if strg[8:].startswith('meanelemp'):
                    dicttemp['indxfixpmeanelem'].append(valu)
            dicttemp['indxfixpmeanelem'] = array(dicttemp['indxfixpmeanelem'])
            dicttemp['indxfixpnumbelemtotl'] = concatenate(dicttemp['indxfixpnumbelem'])
            
            dicttemp['indxfixpdist'] = []
            for strg, valu in dicttemp.iteritems():
                if strg[12:16] == 'dist' and isscalar(valu):
                    dicttemp['indxfixpdist'].append(valu)
                
            dicttemp['indxfixpdist'] = array(dicttemp['indxfixpdist']) 
            dicttemp['indxfixphypr'] = array(list(dicttemp['indxfixpdist']) + list(dicttemp['indxfixpmeanelem']))
        
        if psfnevaltype != 'none':
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    dicttemp['indxfixpsigcene%devt%d' % (i, m)] = cntr.incr()
                    if psfntype == 'doubking' or psfntype == 'singking':
                        dicttemp['indxfixpgamcene%devt%d' % (i, m)] = cntr.incr()
                        if psfntype == 'doubking':
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
        
        if psfnevaltype != 'none' and oaxitype:
            indxfixppsfponor = indxfixppsfpinit + indxpsfponor
            indxfixppsfpoind = indxfixppsfpinit + indxpsfpoind
            indxfixppsfpoaxi = sort(concatenate((indxfixppsfponor, indxfixppsfpoind)))

        dicttemp['indxfixpbacp'] = []
        for c in indxback:
            if specback[c]:
                for d in gdat.indxregi:
                    indx = cntr.incr()
                    dicttemp['indxfixpbacpback%04dreg%d' % (c, d)] = indx
                    dicttemp['indxfixpbacp'].append(indx)
            else:
                for d in gdat.indxregi:
                    for i in gdat.indxener:
                        indx = cntr.incr()
                        dicttemp['indxfixpbacpback%04dreg%dene%d' % (c, d, i)] = indx
                        dicttemp['indxfixpbacp'].append(indx)
                    
        dicttemp['indxfixpbacp'] = array(dicttemp['indxfixpbacp'])
        dicttemp['indxfixphost'] = []
        dicttemp['indxfixpsour'] = []
        dicttemp['indxfixplenp'] = []
        
        # temp
        #dicttemp['indxfixpanglsour'] = []
        #dicttemp['indxfixpanglhost'] = []
        #dicttemp['indxfixpangllens'] = []
        
        if hostemistype != 'none':
            dicttemp['indxfixpspecsour'] = []
            dicttemp['indxfixpspechost'] = []

        for d in gdat.indxregi:
            if lensmodltype != 'none':
                dicttemp['indxfixplgalsourreg%d' % d] = cntr.incr()
                dicttemp['indxfixpbgalsourreg%d' % d] = cntr.incr()
                dicttemp['indxfixpfluxsourreg%d' % d] = cntr.incr()
                if gdat.numbener > 1:
                    dicttemp['indxfixpsindsourreg%d' % d] = cntr.incr()
                dicttemp['indxfixpsizesourreg%d' % d] = cntr.incr()
                dicttemp['indxfixpellpsourreg%d' % d] = cntr.incr()
                dicttemp['indxfixpanglsourreg%d' % d] = cntr.incr()
            if hostemistype != 'none':
                dicttemp['indxfixplgalhostreg%d' % d] = cntr.incr()
                dicttemp['indxfixpbgalhostreg%d' % d] = cntr.incr()
                dicttemp['indxfixpfluxhostreg%d' % d] = cntr.incr()
                if gdat.numbener > 1:
                    dicttemp['indxfixpsindhostreg%d' % d] = cntr.incr()
                dicttemp['indxfixpsizehostreg%d' % d] = cntr.incr()
            if lensmodltype != 'none':
                dicttemp['indxfixpbeinhostreg%d' % d] = cntr.incr()
            if hostemistype != 'none':
                dicttemp['indxfixpellphostreg%d' % d] = cntr.incr()
                dicttemp['indxfixpanglhostreg%d' % d] = cntr.incr()
                dicttemp['indxfixpserihostreg%d' % d] = cntr.incr()
            if lensmodltype != 'none':
                dicttemp['indxfixpsherextrreg%d' % d] = cntr.incr()
                dicttemp['indxfixpsangextrreg%d' % d] = cntr.incr()
                dicttemp['indxfixpsour'] = []
        
        # construct index arrays for individual lens parameters that contain parameters for all regions
        for name, valu in deepcopy(dicttemp).iteritems():
            if name[:-1].endswith('reg'):
                if not name[:-4] in dicttemp:
                    dicttemp[name[:-4]] = []
                dicttemp[name[:-4]].append(valu)
            

        if lensmodltype != 'none' and hostemistype == 'none':
            raise Exception('Lensing cannot be modeled without host galaxy emission.')

        if lensmodltype != 'none':
            liststrgsour = ['lgalsour', 'bgalsour', 'fluxsour', 'sizesour', 'ellpsour', 'anglsour']
            if gdat.numbener > 1:
                liststrgsour += ['sindsour']
        else:
            liststrgsour = []
        
        if hostemistype != 'none':
            liststrghost = ['lgalhost', 'bgalhost', 'fluxhost', 'sizehost', 'ellphost', 'anglhost', 'serihost']
            if gdat.numbener > 1:
                liststrghost += ['sindhost']
            if lensmodltype != 'none':
                liststrghost += ['beinhost', 'sherextr', 'sangextr']

            for strg, valu in dicttemp.iteritems():
                
                if isinstance(valu, list) or isinstance(valu, ndarray):
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
                    
                for d in gdat.indxregi:
                    if strg[8:].startswith('fluxsour') or strg[8:].startswith('sindsour'):
                        dicttemp['indxfixpspecsour'].append(valu)

                    if strg[8:].startswith('fluxhost') or strg[8:].startswith('sindhost'):
                        dicttemp['indxfixpspechost'].append(valu)
                
                if valu in dicttemp['indxfixpsour'] or valu in dicttemp['indxfixphost']:
                    dicttemp['indxfixplenp'].append(valu)
        
        # number of fixed-dimension parameters
        numbfixp = cntr.incr(0)
        
        # indices of fixed-dimension parameters
        indxfixp = arange(numbfixp)
        
        ## number of model spectral parameters for each population
        #numbspep = empty(numbpopl, dtype=int)
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
        indxsamptrap = arange(numbfixp, numbpara)
        if numbtrap > 0:
            indxsampcompinit = indxsamptrap[0]
        else:
            indxsampcompinit = -1
        indxpara = arange(numbpara)

    for strg, valu in locals().iteritems():
        if strg != 'gdat' and '__' not in strg and not strg.endswith('temp') and strg != 'cntr' and len(strg) % 4 == 0:
            setattr(gdat, strgmodl + strg, valu)

    for strg, valu in dicttemp.iteritems():
        # temp
        if isinstance(valu, ndarray) and strg[12:16] != 'dist':
            valu = sort(valu)
        setattr(gdat, strgmodl + strg, valu)


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
    
    elemtype = getattr(gdat, strgmodl + 'elemtype')
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    numbtrapregipoplcuml = getattr(gdat, strgmodl + 'numbtrapregipoplcuml')
    numbtrapregipoplcumr = getattr(gdat, strgmodl + 'numbtrapregipoplcumr')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    numbfixp = getattr(gdat, strgmodl + 'numbfixp')
    indxregipopl = getattr(gdat, strgmodl + 'indxregipopl')
    numbback = getattr(gdat, strgmodl + 'numbback')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    numbpara = getattr(gdat, strgmodl + 'numbpara')
    listnameback = getattr(gdat, strgmodl + 'listnameback')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    lablelemextn = getattr(gdat, strgmodl + 'lablelemextn')
    
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
        if 'pop' in strg:
            
            if numbpopl == 1:
                strgpopl = ''
                strgpoplcomm = ''
            else:
                strgpopl = '%s' % strg[-1]
                strgpoplcomm = ',%s' % strg[-1]
            indxpopltemp = int(strg[-1])
            if namefixp[k].startswith('numbelem'):
                lablfixp[k] = '$N_{%s%s}$' % (lablelemextn[indxpopltemp], strgpoplcomm)
                scalfixp[k] = 'pois'

            if namefixp[k].startswith('meanelem'):
                lablfixp[k] = r'$\mu_{%s%s}$' % (lablelemextn[indxpopltemp], strgpoplcomm)
                scalfixp[k] = scalmeanelem
    
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
                scalfixp[k] = getattr(gdat, strgmodl + 'scal' + namefixp[k])

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
                n = k - getattr(gdat, strgmodl + 'indxfixpsigcene0evt0')
                meanfixp[k] = gdat.psfpexpr[n]
                stdvfixp[k] = 0.1 * gdat.psfpexpr[n]
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
            
            # find the sample vector index of the first background parameter
            try:
                indxfixpbacpinit
            except:
                indxfixpbacpinit = k
            
            # indices of the background and region
            indxbacktemp = int(strgvarb[8:12])
            indxregitemp = int(strgvarb[15])
            
            name = 'bacpback%04dreg%d' % (indxbacktemp, indxregitemp)
            
            if specback[indxbacktemp]:
                strgenertemp = ''
            else:
                i = indxenerbacp[k-indxfixpbacpinit]
                if gdat.numbener > 1:
                    strgenertemp = '%d' % i
                else:
                    strgenertemp = ''
                name += 'ene%d' % i
            
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
                if strgvarb[:-1].endswith('reg'):
                    indxregi = int(strgvarb[-1])
                if gdat.numbregi > 0:
                    strgregi = ',%d' % indxregi
                else:
                    strgregi = ''
                if strgvarb.startswith('lgalsour'):
                    lablfixp[k] = r'$\theta_{\rm{1,src%s}}$' % strgregi
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('bgalsour'):
                    lablfixp[k] = r'$\theta_{\rm{2,src%s}}$' % strgregi
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('fluxsour'):
                    lablfixp[k] = r'$f_{\rm{src}%s}$' % strgregi
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if gdat.numbener > 1:
                    if strgvarb.startswith('sindsour'):
                        lablfixp[k] = r'$s_{\rm{src}%s}$' % strgregi
                        scalfixp[k] = 'self'
                if strgvarb.startswith('sizesour'):
                    lablfixp[k] = r'$\theta_{\rm{e,src}%s}$' % strgregi
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('ellpsour'):
                    lablfixp[k] = r'$\epsilon_{\rm{src}%s}$' % strgregi
                    scalfixp[k] = 'self'
                if strgvarb.startswith('anglsour'):
                    lablfixp[k] = r'$\phi_{\rm{src}%s}$' % strgregi
                    scalfixp[k] = 'self'
                if strgvarb.startswith('lgalhost'):
                    lablfixp[k] = r'$\theta_{1,\rm{hst}%s}$' % strgregi
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('bgalhost'):
                    lablfixp[k] = r'$\theta_{2,\rm{hst}%s}$' % strgregi
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('fluxhost'):
                    lablfixp[k] = r'$f_{\rm{hst}%s}$' % strgregi
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if gdat.numbener > 1:
                    if strgvarb.startswith('sindhost'):
                        lablfixp[k] = r'$s_{\rm{hst}%s}$' % strgregi
                        scalfixp[k] = 'self'
                if strgvarb.startswith('sizehost'):
                    lablfixp[k] = r'$\theta_{\rm{e,hst}%s}$' % strgregi
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('beinhost'):
                    lablfixp[k] = r'$\theta_{\rm{E,hst}%s}$' % strgregi
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb.startswith('ellphost'):
                    lablfixp[k] = r'$\epsilon_{\rm{hst}%s}$' % strgregi
                    scalfixp[k] = 'self'
                if strgvarb.startswith('anglhost'):
                    lablfixp[k] = r'$\phi_{\rm{hst}%s}$' % strgregi
                    scalfixp[k] = 'self'
                if strgvarb.startswith('serihost'):
                    lablfixp[k] = r'$n_{\rm{S,hst}%s}$' % strgregi
                    scalfixp[k] = 'self'
                if strgvarb.startswith('sherextr'):
                    lablfixp[k] = r'$\gamma_{ext%s}$' % strgregi
                    scalfixp[k] = 'self'
                if strgvarb.startswith('sangextr'):
                    lablfixp[k] = r'$\phi_{ext%s}$' % strgregi
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
            factfixp[k] = log(maxmfixp[k] / minmfixp[k])
        if scalfixp[k] == 'atan':
            factfixp[k] = arctan(maxmfixp[k]) - arctan(minmfixp[k])
        if scalfixp[k] == 'gaus':
            minmfixp[k] = meanfixp[k] - 3. * stdvfixp[k]
            maxmfixp[k] = meanfixp[k] + 3. * stdvfixp[k]
        if scalfixp[k] == 'eerr':
            cdfnminmfixp[k], cdfndifffixp[k] = retr_eerrnorm(minmfixp[k], maxmfixp[k], meanfixp[k], stdvfixp[k])
   
    # background templates
    listlegdsbrt = deepcopy(listlegdback)
    if 'lght' in elemtype:
        listlegdsbrt.append('PS')
        listlegdsbrt.append('PS Unres')
    if 'lens' in elemtype:
        listlegdsbrt.append('Source')
        listlegdsbrt.append('Host')
    if 'clus' in elemtype:
        listlegdsbrt.append('Uniform')
    listlegdsbrtspec = ['Data']
    listlegdsbrtspec += deepcopy(listlegdsbrt)
    if len(listlegdsbrt) > 1:
        listlegdsbrtspec.append('Total Model')
    numblablsbrt = len(listlegdsbrt)
    numblablsbrtspec = len(listlegdsbrtspec)
    
    namepara = zeros(numbpara, dtype=object)
    lablpara = zeros(numbpara, dtype=object)
    scalpara = zeros(numbpara, dtype=object)
    factplotpara = zeros(numbpara)
    
    listfactplotcomp = []
    listlablcomp = []
    for l in indxpopl:
        listfactplotcomp.append(zeros(numbcomp[l]))
        listlablcomp.append(zeros(numbcomp[l], dtype=object))
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
            for d in indxregipopl[l]:
                if k >= numbtrapregipoplcuml[l][d]:
                    indxpopltemp = l
                    indxregitemp = d
                    indxelemtemp = (k - numbtrapregipoplcuml[indxpopltemp][indxregitemp]) // numbcomp[indxpopltemp]
                    indxcomptemp = (k - numbtrapregipoplcuml[indxpopltemp][indxregitemp]) % numbcomp[indxpopltemp]
                    break
        namepara[numbfixp+k] = '%spop%dreg%d%04d' % (liststrgcomp[indxpopltemp][indxcomptemp], indxpopltemp, indxregitemp, indxelemtemp)
        lablpara[numbfixp+k] = listlablcomp[indxpopltemp][indxcomptemp]
        scalpara[numbfixp+k] = listscalcomp[indxpopltemp][indxcomptemp]
        factplotpara[numbfixp+k] = listfactplotcomp[indxpopltemp][indxcomptemp]
    
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
            if len(strg) % 4 == 0:
                setattr(gdat, strgmodl + strg, valu)


def setp_namevarbsing(gdat, strgmodl, strgvarb, popl, ener, evtt, back, regi):
    
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
    elif back != 'none':
        indxbacktemp = [back]
    
    if regi == 'full':
        indxregitemp = gdat.indxregi
    elif regi != 'none':
        indxregitemp = [regi]
    
    liststrgvarb = []
    if popl != 'none' and ener == 'none' and evtt == 'none' and back == 'none' and regi == 'none':
        for l in indxpopltemp:
            liststrgvarb.append(strgvarb + 'pop%d' % l)
    if popl != 'none' and ener == 'none' and evtt == 'none' and back == 'none' and regi != 'none':
        for d in indxregitemp:
            for l in indxpopltemp:
                liststrgvarb.append(strgvarb + 'pop%dreg%d' % (l, d))
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none' and regi != 'none':
        for d in indxregitemp:
            liststrgvarb.append(strgvarb + 'reg%d' % d)
    if popl == 'none' and ener != 'none' and evtt != 'none' and back == 'none' and regi == 'none':
        for i in indxenertemp:
            for m in indxevtttemp:
                liststrgvarb.append(strgvarb + 'ene%devt%d' % (i, m))
    if popl == 'none' and ener != 'none' and evtt == 'none' and back != 'none' and regi != 'none':
        for c in indxbacktemp:
            for d in indxregitemp:
                for i in indxenertemp:
                    liststrgvarb.append(strgvarb + 'back%04dreg%dene%d' % (c, d, i))
    if popl == 'none' and ener == 'none' and evtt == 'none' and back != 'none' and regi != 'none':
        for c in indxbacktemp:
            for d in indxregitemp:
                liststrgvarb.append(strgvarb + 'back%04dreg%d' % (c, d))
    if popl == 'none' and ener != 'none' and evtt == 'none' and back != 'none' and regi == 'none':
        for c in indxbacktemp:
            for i in indxenertemp:
                liststrgvarb.append(strgvarb + 'back%04dene%d' % (c, i))
    if popl == 'none' and ener != 'none' and evtt == 'none' and back == 'none' and regi == 'none':
        for i in indxenertemp:
            liststrgvarb.append(strgvarb + 'ene%d' % i)
    if popl == 'none' and ener == 'none' and evtt == 'none' and back == 'none' and regi == 'none':
        liststrgvarb.append(strgvarb)
    
    return liststrgvarb


def setp_varbvalu(gdat, strgvarb, valu, popl='none', ener='none', evtt='none', back='none', regi='none', strgmodl=None):
    
    if strgmodl == None:
        if gdat.datatype == 'mock':
            liststrgmodl = ['true', 'fitt']
        else:
            liststrgmodl = ['fitt']
    else:
        liststrgmodl = [strgmodl]
    liststrgmodl = deepcopy(liststrgmodl)
    for strgmodltemp in liststrgmodl:
        liststrgvarb = setp_namevarbsing(gdat, strgmodltemp, strgvarb, popl, ener, evtt, back, regi)
        for strgvarbtemp in liststrgvarb:
            setp_varbiter(gdat, strgmodltemp, strgvarbtemp, valu)


def setp_varbiter(gdat, strgmodltemp, strgvarbtemp, valu):
    try:
        valutemp = getattr(gdat, strgvarbtemp)
        if valutemp == None:
            raise
        setattr(gdat, strgmodltemp + strgvarbtemp, valutemp)
    except:
        try:
            valutemp = getattr(gdat, strgmodltemp + strgvarbtemp)
            if valutemp == None:
                raise
        except:
            setattr(gdat, strgmodltemp + strgvarbtemp, valu)
 

def setp_varblimt(gdat, strgvarb, listvalu, typelimt='minmmaxm', popl='none', ener='none', evtt='none', back='none', regi='none', strgmodl=None):
    
    if strgmodl == None:
        if gdat.datatype == 'mock':
            liststrgmodl = ['true', 'fitt']
        else:
            liststrgmodl = ['fitt']
    else:
        liststrgmodl = [strgmodl]
    
    for strgmodltemp in liststrgmodl:

        liststrgvarb = setp_namevarbsing(gdat, strgmodltemp, strgvarb, popl, ener, evtt, back, regi)

        for strgvarbtemp in liststrgvarb:

            if typelimt == 'minmmaxm':
                setp_varbiter(gdat, strgmodltemp, 'minm' + strgvarbtemp, listvalu[0])
                setp_varbiter(gdat, strgmodltemp, 'maxm' + strgvarbtemp, listvalu[1])
            else:
                setp_varbiter(gdat, strgmodltemp, 'mean' + strgvarbtemp, listvalu[0])
                setp_varbiter(gdat, strgmodltemp, 'stdv' + strgvarbtemp, listvalu[1])
                
                # set minimum and maximum for Gaussian distributed variables
                meanpara = getattr(gdat, strgmodltemp + 'mean' + strgvarbtemp)
                stdvpara = getattr(gdat, strgmodltemp + 'stdv' + strgvarbtemp)
                setattr(gdat, strgmodltemp + 'minm' + strgvarbtemp, meanpara - stdvpara * 5)
                setattr(gdat, strgmodltemp + 'maxm' + strgvarbtemp, meanpara + stdvpara * 5)
 

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


def init_figr(gdat, gdatmodi, strgplot, strgstat, strgmodl, indxregiplot, indxenerplot, indxevttplot, indxpoplplot, intreval=False):
    
    if intreval:
        figrsize = [10, 10]
    else:
        figrsize = (gdat.sizeimag, gdat.sizeimag)
    figr, axis = plt.subplots(figsize=figrsize)
    
    if intreval:
        axis.set_position([0.5, 0.1, 0.6, 0.8])
        path = ''
    else:
        strgregi = 'reg%d' % indxregiplot
        
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

        nameplot = '%s%s%s%s%s' % (strgplot, strgregi, strgener, strgevtt, strgpopl)
   
        path = retr_plotpath(gdat, gdatmodi, strgstat, strgmodl, nameplot)
    
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
    axis.axvline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axvline(-innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(-innr, ls='--', alpha=gdat.alphbndr, color='black')


def retr_imag(gdat, axis, maps, strgstat, strgmodl, strgcbar, indxregiplot=None, indxenerplot=None, indxevttplot=-1, tdim=False, imag=None):
    
    vmin = getattr(gdat, 'minmscal' + strgcbar)
    vmax = getattr(gdat, 'maxmscal' + strgcbar)
    cmap = getattr(gdat, 'cmap' + strgcbar) 
    scal = getattr(gdat, 'scal' + strgcbar) 

    draw_frambndr(gdat, axis)
    
    # flatten the array
    #if tdim:
    #    if indxenerplot == None:
    #        maps = maps.flatten()
    #    else:
    #        maps = maps.reshape((gdat.numbregi, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    
    # take the relevant energy and PSF bins
    if indxenerplot != None:
        if indxevttplot == -1:
            maps = sum(maps[indxregiplot][indxenerplot, ...], axis=1)
        else:
            maps = maps[indxregiplot][indxenerplot, :, indxevttplot]
    else:
        maps = maps[indxregiplot, :]
    
    # project the map to 2D
    if gdat.pixltype == 'heal':# and not tdim:
        maps = tdpy.util.retr_cart(maps, indxpixlrofi=gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                                                                            minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                            minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata)
    
    if gdat.pixltype == 'cart':# or tdim:
        shap = [gdat.numbsidecart] + list(maps.shape)
        shap[1] = gdat.numbsidecart
        shapflat = list(maps.shape)
        shapflat[0] = gdat.numbpixlfull
        mapstemp = zeros(shapflat)
        mapstemp[gdat.indxpixlrofi, ...] = maps
        maps = mapstemp.reshape(shap).swapaxes(0, 1)
    
    # rescale the map
    if scal == 'asnh':
        maps = arcsinh(maps)
    if scal == 'logt':
        maps = log10(maps)
    if imag == None:
        imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest', vmin=vmin, vmax=vmax, alpha=gdat.alphmaps)
        return imag
    else:
        imag.set_data(maps)
    

def make_cbar(gdat, axis, imag, indxenerplot=None, tick=None, labl=None):

    # make a color bar
    cbar = plt.colorbar(imag, ax=axis, fraction=0.05, aspect=15)
    if tick != None and labl != None:
        cbar.set_ticks(tick)
        cbar.set_ticklabels(labl)
    
    return cbar


def make_catllabl(gdat, strgstat, strgmodl, axis):
    
    # transdimensional elements
    if strgmodl == 'fitt' and (strgstat == 'post' and gdat.condcatl or strgstat == 'this') and gdat.fittnumbtrap > 0:
        for l in gdat.fittindxpopl:
            if strgstat == 'post':
                colr = 'black'
                labl = 'Condensed Model %s' % gdat.fittlistlegdelem[l]
            else:
                colr = 'b'
                labl = 'Sample Model %s' % gdat.fittlistlegdelem[l]
            if not gdat.fittmaxmnumbelempopl[l] == 0:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                                                                                label=labl, marker='+', lw=gdat.mrkrlinewdth, color=colr)
    
    if gdat.allwrefr:
        for q in gdat.indxrefr:
            if not amax(gdat.refrnumbelem[q]) == 0:
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                                                                      label=gdat.legdrefrhits[q], marker='x', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
                axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                      label=gdat.legdrefrmiss[q], marker='s', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
    
    # fixed-dimensional objects
    if strgmodl == 'fitt':
        if gdat.fittlensmodltype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='Model Source', marker='<', lw=gdat.mrkrlinewdth, color='b')
        
        if gdat.fitthostemistype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='Model Host', marker='s', lw=gdat.mrkrlinewdth, color='b')
    if strgmodl == 'true':
        if gdat.truelensmodltype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='%s Source' % gdat.legdrefr[q], marker='>', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
        
        if gdat.truehostemistype != 'none':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='%s Host' % gdat.legdrefr[q], marker='D', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
    
    temphand, temp = axis.get_legend_handles_labels()
    numblabl = len(temp)
    
    if numblabl == 4:
        numbcols = 2
    else:
        numbcols = 3
    axis.legend(bbox_to_anchor=[0.5, 1.15], loc='center', ncol=numbcols)
        

def supr_fram(gdat, gdatmodi, strgstat, strgmodl, axis, indxregiplot, indxpoplplot=-1):
    
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    
    # associations with the reference elements
    if gdat.allwrefr:
        for q in gdat.indxrefr:
            if gdat.refrnumbelem[q][indxregiplot] > 0:
                if indxpoplplot == -1:
                    listindxpoplplot = gdat.fittindxpopl
                else:
                    listindxpoplplot = [indxpoplplot]
                for l in listindxpoplplot:
                    reframpl = getattr(gdat, 'refr' + gdat.listnamefeatamplrefr[q])
                    if reframpl == None:
                        mrkrsize = full(gdat.refrnumbelem[q][indxregiplot], 5.)
                    else:
                        mrkrsize = retr_mrkrsize(gdat, reframpl[q][indxregiplot][0, :], gdat.listnamefeatamplrefr[q])
                    lgal = copy(gdat.refrlgal[q][indxregiplot][0, :])
                    bgal = copy(gdat.refrbgal[q][indxregiplot][0, :])
                    numbelem = int(gdat.refrnumbelem[q][indxregiplot])

                    if gdatmodi != None and numbtrap > 0:   
                        ### hit
                        indx = gdatmodi.thisindxelemrefrasschits[q][l][indxregiplot]
                        if indx.size > 0:
                            axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, label=gdat.legdrefrmiss, \
                                                                                                                      marker='x', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
                        ### missed
                        indx = gdatmodi.thisindxelemrefrasscmiss[q][l][indxregiplot]
                    else:
                        indx = arange(lgal.size)
                    if indx.size > 0: 
                        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, facecolor='none', \
                                                                                            label=gdat.legdrefrhits, marker='s', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
            
        # temp -- generalize this to input refrlgalhost vs.
        if gdat.datatype == 'mock':
            ## host galaxy position
            if hostemistype != 'none':
                axis.scatter(gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixplgalhost[indxregiplot]], \
                                                                            gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixpbgalhost[indxregiplot]], \
                                                                            facecolor='none', \
                                                                            #alpha=gdat.alphelem, \
                                                                            alpha=0.7, \
                                                                            label=gdat.legdrefrhits, s=300, marker='D', lw=gdat.mrkrlinewdth, color='g')
            if lensmodltype != 'none':
                ## host galaxy Einstein radius
                axis.add_patch(plt.Circle((gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixplgalhost[indxregiplot]], \
                                                                                gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixpbgalhost[indxregiplot]]), \
                                                                                gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixpbeinhost[indxregiplot]], \
                                                                                edgecolor='g', facecolor='none', lw=gdat.mrkrlinewdth))
                
                ## source galaxy position
                axis.scatter(gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixplgalsour[indxregiplot]], \
                                                                            gdat.anglfact * gdat.truesampvarb[gdat.trueindxfixpbgalsour[indxregiplot]], \
                                                                            facecolor='none', \
                                                                            alpha=0.7, \
                                                                            #alpha=gdat.alphelem, \
                                                                            label=gdat.legdrefrhits, s=300, marker='>', lw=gdat.mrkrlinewdth, color='g')
        
    # model catalog
    if indxpoplplot == -1:
        listindxpoplplot = gdat.fittindxpopl
    else:
        listindxpoplplot = [indxpoplplot]
    for l in listindxpoplplot:
        if gdatmodi != None:
            if gdat.fittnumbtrap > 0:
                mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.fittnamefeatampl[l]][l][indxregiplot]], gdat.fittnamefeatampl[l])
                lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][l][indxregiplot]]
                bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][l][indxregiplot]]
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphelem, label='Sample', marker='+', lw=gdat.mrkrlinewdth, color='b')

            ## source
            if lensmodltype != 'none':
                lgalsour = gdatmodi.thissampvarb[gdat.fittindxfixplgalsour[indxregiplot]]
                bgalsour = gdatmodi.thissampvarb[gdat.fittindxfixpbgalsour[indxregiplot]]
                axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, facecolor='none', \
                                                                      alpha=gdat.alphelem, \
                                                                      label='Model Source', s=300, marker='<', lw=gdat.mrkrlinewdth, color='b')
    
            if hostemistype != 'none':
                ## host
                lgalhost = gdatmodi.thissampvarb[gdat.fittindxfixplgalhost[indxregiplot]]
                bgalhost = gdatmodi.thissampvarb[gdat.fittindxfixpbgalhost[indxregiplot]]
                axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, facecolor='none', \
                                                                      alpha=gdat.alphelem, \
                                                                      label='Model Host', s=300, marker='s', lw=gdat.mrkrlinewdth, color='b')
                if lensmodltype != 'none':
                    beinhost = gdatmodi.thissampvarb[gdat.fittindxfixpbeinhost[indxregiplot]]
                    axis.add_patch(plt.Circle((gdat.anglfact * lgalhost, gdat.anglfact * bgalhost), gdat.anglfact * beinhost, edgecolor='b', facecolor='none', \
                                                                                                                                         lw=gdat.mrkrlinewdth, ls='--'))
                
    # temp
    if strgstat == 'post' and gdat.condcatl and gdat.fittnumbtrap > 0:
        lgal = zeros(gdat.numbprvlhigh)
        bgal = zeros(gdat.numbprvlhigh)
        ampl = zeros(gdat.numbprvlhigh)
        cntr = 0
        for r in gdat.indxstkscond:
            if r in gdat.indxprvlhigh:
                lgal[cntr] = gdat.dictglob['poststkscond'][r]['lgal'][0]
                bgal[cntr] = gdat.dictglob['poststkscond'][r]['bgal'][0]
                # temp -- this does not allow sources with different spectra to be assigned to the same stacked sample
                ampl[cntr] = gdat.dictglob['poststkscond'][r][gdat.fittnamefeatampl[l]][0]
                cntr += 1
        mrkrsize = retr_mrkrsize(gdat, ampl, gdat.fittnamefeatampl[l])
        axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, label='Condensed', marker='+', color='black', lw=gdat.mrkrlinewdth)
        if False:
            for r in gdat.indxstkscond:
                lgal = array([gdat.dictglob['liststkscond'][r]['lgal']])
                bgal = array([gdat.dictglob['liststkscond'][r]['bgal']])
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, marker='+', color='black', alpha=0.1, lw=gdat.mrkrlinewdth)


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


def writfile(gdattemp, path):
    
    filepick = open(path + '.p', 'wb')
    filearry = h5py.File(path + '.h5', 'w')
    
    gdattemptemp = tdpy.util.gdatstrt()
    for attr, valu in gdattemp.__dict__.iteritems():
        
        if isinstance(valu, ndarray) and valu.dtype != dtype('O') or isinstance(valu, str) or \
                                                        isinstance(valu, float) or isinstance(valu, bool) or isinstance(valu, int) or isinstance(valu, float64):
            filearry.create_dataset(attr, data=valu)
        else:
            setattr(gdattemptemp, attr, valu)

    cPickle.dump(gdattemptemp, filepick, protocol=cPickle.HIGHEST_PROTOCOL)
    filepick.close()
    filearry.close()
   

def readoutp(path):
    
    thisfile = h5py.File(path, 'r')
    gdattemp = tdpy.util.gdatstrt()
    for attr in thisfile:
        setattr(gdattemp, attr, thisfile[attr]) 
    thisfile.close()

    return gdattemp


def writoutp(gdat, path):

    thisfile = h5py.File(path, 'w')
    for attr, valu in gdat.__dict__.iteritems():
        if attr.startswith('list'):
            attrtemp = attr[4:]
            if attrtemp in gdat.fittliststrgfeatodimtotl:
                for l in gdat.fittindxpopl:
                    for d in gdat.fittindxregipopl[l]:
                        attrprim = attrtemp + 'pop%dreg%d' % (l, d)
                        listtemp = []
                        arrytemp = zeros((gdat.numbsamptotl, gdat.fittmaxmnumbelem[l][d])) + nan
                        for n in range(len(valu)):
                            arrytemp[n, :len(valu[n][l][d])] = valu[n][l][d]
                        thisfile.create_dataset(attrprim, data=arrytemp)
            elif attrtemp == 'meanelem' or attrtemp == 'fluxdistslop' or attrtemp == 'psfp' or attrtemp == 'bacp':
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


def initchro(gdat, gdatmodi, strgstat, name):

    if strgstat == 'next':
        gdatmodi.thischro[gdat.indxchro[name]] = gdat.functime()
    if strgstat == 'gene':
        gdatmodi.thischro[name] = gdat.functime()
    

def stopchro(gdat, gdatmodi, strgstat, name):
    
    if strgstat == 'next':
        gdatmodi.thischro[gdat.indxchro[name]] = gdat.functime() - gdatmodi.thischro[gdat.indxchro[name]]
    if strgstat == 'gene':
        gdatmodi.thischro[name] = gdat.functime() - gdatmodi.thischro[name]


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


def retr_lpriselfdist(gdat, strgmodl, feat, strgfeat):
    
    minm = getattr(gdat, 'minm' + strgfeat)
    maxm = getattr(gdat, 'maxm' + strgfeat)
    
    lpri = sum(log(pdfn_self(feat, minm, maxm)))
    
    return lpri


def retr_lpripowrdist(gdat, gdatmodi, strgmodl, feat, strgfeat, sampvarb, l):
    
    minm = getattr(gdat, strgmodl + 'minm' + strgfeat)
    maxm = getattr(gdat, strgmodl + 'maxm' + strgfeat)
    distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
    
    lpri = sum(log(pdfn_powr(feat, minm, maxm, distslop)))
    
    return lpri


def retr_lprigausdist(gdat, gdatmodi, strgmodl, feat, strgfeat, sampvarb, l):
    
    distmean = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distmean')[l]]
    diststdv = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'diststdv')[l]]
    
    lpri = sum(log(pdfn_gaus(feat, distmean, diststdv)))
    
    return lpri


def retr_lpriigamdist(gdat, gdatmodi, strgmodl, feat, strgfeat, sampvarb, l):
    
    distslop = sampvarb[getattr(gdat, strgmodl + 'indxfixp' + strgfeat + 'distslop')[l]]
    cutf = getattr(gdat, strgmodl + 'cutf' + strgfeat)
    
    lpri = sum(log(pdfn_igam(feat, distslop, cutf)))

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
def retr_deflelem_jitt(deflelem, indxpixl, lgalgrid, bgalgrid, numbelemconc, lgalconc, bgalconc, defsconc, ascaconc, acutconc):
    
    for k in range(numbelemconc):
        deflelem[:] += retr_defl_jitt(indxpixl, lgalgrid, bgalgrid, lgalconc[k], bgalconc[k], defsconc[k], 0., 0., asca=ascaconc[k], acut=acutconc[k])


def retr_strgpfix(strgstat, strgmodl):
    
    if strgmodl == 'fitt':
        strgpfix = strgstat
    if strgmodl == 'true':
        if strgstat == 'this':
            strgpfix = 'true'
        else:
            strgpfix = strgstat + strgmodl

    return strgpfix

        
def retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl):
    
    if strgmodl == 'true' or strgmodl == 'fitt' and (strgstat == 'mlik' or strgstat == 'post'):
        gdatobjt = gdat
    else:
        gdatobjt = gdatmodi

    return gdatobjt


def proc_samp(gdat, gdatmodi, strgstat, strgmodl, raww=False, fast=False):
    
    initchro(gdat, gdatmodi, strgstat, 'proc')

    if gdat.verbtype > 1:
        print 'proc_samp()'
        print 'strgstat'
        print strgstat
        print
    
    gdatobjt = retr_gdatobjt(gdat, gdatmodi, strgstat, strgmodl)
    
    fittnumbelemzero = getattr(gdat, strgmodl + 'numbelemzero')
    # temp
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        indxfixpmeanelem = getattr(gdat, strgmodl + 'indxfixpmeanelem')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        elemtype = getattr(gdat, strgmodl + 'elemtype')
        boolelemspectotl = getattr(gdat, strgmodl + 'boolelemspectotl')
        if boolelemspectotl:
            minmflux = getattr(gdat, strgmodl + 'minmflux')
        if 'lens' in elemtype:
            minmdefs = getattr(gdat, strgmodl + 'minmdefs')
        if 'clus' in elemtype:
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
    
        if gdat.verbtype > 1:
            print 'spectype'
            print spectype
    
    backtype = getattr(gdat, strgmodl + 'backtype')
    psfntype = getattr(gdat, strgmodl + 'psfntype')
    convdiff = getattr(gdat, strgmodl + 'convdiff')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    listnamediff = getattr(gdat, strgmodl + 'listnamediff')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')
    
    strgpfix = retr_strgpfix(strgstat, strgmodl)
    strgpfixthis = retr_strgpfix('this', strgmodl)

    # grab the sample vector
    sampvarb = getattr(gdatobjt, strgpfix + 'sampvarb')
    
    if gdat.diagmode:
        if not isfinite(sampvarb).all():
            raise Exception('')

    if psfnevaltype != 'none':
        psfp = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfp')]
        oaxitype = getattr(gdat, strgmodl + 'oaxitype')
    
    bacp = sampvarb[getattr(gdat, strgmodl + 'indxfixpbacp')]
    
    if numbtrap > 0:
        indxregipopl = getattr(gdat, strgmodl + 'indxregipopl')
        indxelemfull = list(getattr(gdatobjt, strgpfix + 'indxelemfull'))
        # temp -- this may slow down execution
        indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, strgmodl)
        setattr(gdatobjt, strgpfix + 'indxsampcomp', indxsampcomp)
        indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
    
        numbelempopl = empty(numbpopl)
        numbelem = [[] for l in indxpopl]
        for l in indxpopl:
            numbelem[l] = sampvarb[indxfixpnumbelem[l]].astype(int)
            numbelempopl[l] = sum(numbelem[l])
        numbelemtotl = sum(numbelempopl) 

        dictelem = [[[] for d in indxregipopl[l]] for l in indxpopl]
        for l in indxpopl:
            for d in indxregipopl[l]:
                dictelem[l][d] = dict()
                for strgfeat in liststrgfeatdefa:
                    dictelem[l][d][strgfeat] = []
                for strgcomp in liststrgcomp[l]:
                    dictelem[l][d][strgcomp] = sampvarb[indxsampcomp[strgcomp][l][d]]
        
        if gdat.diagmode:
            for l in indxpopl:
                for d in indxregipopl[l]:
                    if numbelem[l][d] != dictelem[l][d]['lgal'].size:
                        print 'numbelem'
                        print numbelem
                        print 'dictelem[0][0][lgal]'
                        print dictelem[l][d]['lgal']
                        raise Exception('')

        if gdat.verbtype > 1:
            for l in indxpopl:
                print 'l'
                print l
                for d in indxregipopl[l]:
                    print 'd'
                    print d
                    print 'dictelem[l][d]'
                    for strgcomp in liststrgcomp[l]:
                        print strgcomp
                        print dictelem[l][d][strgcomp]
                print
        if gdat.diagmode:
            for l in indxpopl:
                for d in indxregipopl[l]:
                    for strgcomp in liststrgcomp[l]:
                        if (dictelem[l][d][strgcomp] < getattr(gdat, strgmodl + 'minm' + strgcomp)).any() or \
                                            (dictelem[l][d][strgcomp] > getattr(gdat, strgmodl + 'maxm' + strgcomp)).any():
                            raise Exception('Element parameter outside prior.')
           
        # temp
        if False and sum(numbelem) == 1:
            try:
                gdatmodi.cntrswep
                booltemp = True
            except:
                booltemp = False
            if booltemp:
                print 'dictelem[0][0][lgal]'
                print dictelem[0][0]['lgal']
                print 'dictelem[0][0][bgal]'
                print dictelem[0][0]['bgal']
                if abs(dictelem[0][0]['lgal']) / gdat.maxmgang > 1e-9 or abs(dictelem[0][0]['bgal']) / gdat.maxmgang > 1e-9:
                    raise Exception('')
                print

        # temp
        if gdat.diagmode:
            for l in indxpopl:
                if elemtype[l] == 'lens':
                    for d in indxregipopl[l]:
                        if gdat.variasca:
                            indx = where(sampvarb[indxsampcomp['acut'][l][d]] < 0.)[0]
                            if indx.size > 0:
                                print 'Acut went negative'
                                sampvarb[indxsampcomp['acut'][l][d]][indx] = 1e-3 * gdat.anglfact
                                raise Exception('')
                        if gdat.variacut:
                            indx = where(sampvarb[indxsampcomp['asca'][l][d]] < 0.)[0]
                            if indx.size > 0:
                                print 'Asca went negative'
                                sampvarb[indxsampcomp['asca'][l][d]][indx] = 1e-3 * gdat.anglfact
                                raise Exception('')

        if strgstat == 'this' or strgstat == 'next' and (('lght' in elemtype or 'clus' in elemtype) and gdatmodi.proppsfp):
            # evaluate spectra
            for l in indxpopl:
                if elemtype[l] == 'lght' or elemtype[l] == 'lghtline':
                    for d in indxregipopl[l]:
                        if elemtype[l] == 'lght':
                            dictelem[l][d]['spec'] = retr_spec(gdat, dictelem[l][d]['flux'], sind=dictelem[l][d]['sind'], curv=dictelem[l][d]['curv'], \
                                               expc=dictelem[l][d]['expc'], sind0001=dictelem[l][d]['sind0001'], sind0002=dictelem[l][d]['sind0002'], spectype=spectype[l])
                        else:
                            dictelem[l][d]['spec'] = retr_spec(gdat, dictelem[l][d]['flux'], elin=dictelem[l][d]['elin'], edis=gdat.edis, spectype=spectype[l])
            
    ### loglikelihood
    initchro(gdat, gdatmodi, strgstat, 'llik')
    
    # temp -- this neglects prior terms that depend on the derived quantities
    evalllik = False
    if strgstat == 'next': 
        if gdatmodi.propllik:
            indxregieval = gdatmodi.indxregimodi
            indxenereval = gdatmodi.indxenermodi
            indxevtteval = gdatmodi.indxevttmodi
            evalllik = True
    else:
        indxregieval = gdat.indxregi 
        indxenereval = gdat.indxener 
        indxevtteval = gdat.indxevtt
        evalllik = True
    
    if evalllik or strgstat != 'next':
        
        if gdat.verbtype > 1:
            print 'Evaluating the likelihood...'
        
        # process a sample vector and the occupancy list to calculate secondary variables
        if lensmodltype != 'none':
            lgalsour = sampvarb[getattr(gdat, strgmodl + 'indxfixplgalsour')]
            bgalsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpbgalsour')]
            fluxsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpfluxsour')]
            if gdat.numbener > 1:
                sindsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpsindsour')]
            sizesour = sampvarb[getattr(gdat, strgmodl + 'indxfixpsizesour')]
            ellpsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpellpsour')]
            anglsour = sampvarb[getattr(gdat, strgmodl + 'indxfixpanglsour')]
        if hostemistype != 'none':
            lgalhost = sampvarb[getattr(gdat, strgmodl + 'indxfixplgalhost')]
            bgalhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpbgalhost')]
            fluxhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpfluxhost')]
            if gdat.numbener > 1:
                sindhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpsindhost')]
            sizehost = sampvarb[getattr(gdat, strgmodl + 'indxfixpsizehost')]
        if lensmodltype != 'none':
            if raww:
                beinhost = zeros(gdat.numbregi)
            else:
                beinhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpbeinhost')]
        if hostemistype != 'none':
            ellphost = sampvarb[getattr(gdat, strgmodl + 'indxfixpellphost')]
            anglhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpanglhost')]
            serihost = sampvarb[getattr(gdat, strgmodl + 'indxfixpserihost')]
        if lensmodltype != 'none':
            initchro(gdat, gdatmodi, strgstat, 'deflzero')
            defl = zeros((gdat.numbregi, gdat.numbpixl, 2))
            stopchro(gdat, gdatmodi, strgstat, 'deflzero')
            
        if lensmodltype != 'none':
            if raww:
                sherextr = zeros(gdat.numbregi)
            else:
                sherextr = sampvarb[getattr(gdat, strgmodl + 'indxfixpsherextr')]
            sangextr = sampvarb[getattr(gdat, strgmodl + 'indxfixpsangextr')]
           
            ## host halo deflection
            initchro(gdat, gdatmodi, strgstat, 'deflhost')
            if strgstat == 'next' and not gdatmodi.prophost:
                # retrieve state variable
                deflhost = getattr(gdatobjt, strgpfixthis + 'deflhost')
            else:
                deflhost = [[] for d in indxregieval]
                for dd, d in enumerate(indxregieval):
                    deflhost[dd] = retr_defl(gdat, lgalhost[d], bgalhost[d], beinhost[d], ellphost[d], anglhost[d])
                setattr(gdatobjt, strgpfix + 'deflhost', deflhost)
            
            if gdat.diagmode:
                if not isfinite(deflhost).all():
                    raise Exception('')
            
            defl += deflhost

            stopchro(gdat, gdatmodi, strgstat, 'deflhost')

            ## external shear
            initchro(gdat, gdatmodi, strgstat, 'deflextr')
            deflextr = [[] for d in indxregieval]
            for dd, d in enumerate(indxregieval):
                deflextr[dd] = retr_deflextr(gdat, sherextr[d], sangextr[d])
            defl += deflextr
            stopchro(gdat, gdatmodi, strgstat, 'deflextr')
        
        ## construct the PSF to be convolved with the image
        if psfnevaltype == 'conv' or psfnevaltype == 'full':
            initchro(gdat, gdatmodi, strgstat, 'psfnconv')
            if strgstat == 'next' and not gdatmodi.proppsfnconv:
                # retrieve state variable
                psfnconv = getattr(gdatobjt, strgpfixthis + 'psfnconv')
            else:
                psfnconv = [[[] for i in gdat.indxener] for m in gdat.indxevtt]
                if gdat.pixltype == 'cart':
                    if psfntype != 'singgaus':
                        psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
                        fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
                    for m in gdat.indxevtt:
                        for i in gdat.indxener:
                            if psfntype == 'singgaus':
                                sigm = psfp[i+m*gdat.numbener]
                            else:
                                sigm = fwhm[i, m] / 2.355
                            psfnconv[m][i] = AiryDisk2DKernel(sigm / gdat.sizepixl)
                setattr(gdatobjt, strgpfix + 'psfp', psfp)
                setattr(gdatobjt, strgpfix + 'psfnconv', psfnconv)
            stopchro(gdat, gdatmodi, strgstat, 'psfnconv')
    
        if numbtrap > 0 and ('lght' in elemtype or 'clus' in elemtype):
            ## PSF off-axis factor
            if oaxitype:
                onor = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfponor')]
                oind = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfpoind')]
                factoaxi = retr_factoaxi(gdat, gdat.binsoaxi, onor, oind)
        
            psfntype = getattr(gdat, strgmodl + 'psfntype')
            
            if gdat.kernevaltype == 'ulip':
                psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
                if oaxitype:
                    psfnintp = []
                    for p in gdat.indxoaxi:
                        psfnintp.append(interp1d_pick(gdat.binsangl, psfn[:, :, :, p], axis=1))
                else:
                    psfnintp = interp1d_pick(gdat.binsangl, psfn, axis=1)

            if gdat.kernevaltype == 'bspx':
                
                psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsanglcart.flatten(), psfntype, gdat.binsoaxi, oaxitype)
                
                # side length of the upsampled kernel
                gdat.numbsidekernusam = 100
                # side length of the original kernel
                gdat.numbsidekern = gdat.numbsidekernusam / factkernusam 
                gdat.indxsidekern = arange(gdat.numbsidekern)

        		# pad by one row and one column
        		#psf = zeros((gdat.numbsidekernusam+1, gdat.numbsidekernusam+1))
        		#psf[0:gdat.numbsidekernusam, 0:gdat.numbsidekernusam] = psf0
				
        		# make design matrix for each factkernusam x factkernusam region
                nx = factkernusam + 1
                y, x = mgrid[0:nx, 0:nx] / float(factkernusam)
                x = x.flatten()
                y = y.flatten()
                print 'x'
                print x
                kernmatrdesi = array([full(nx*nx, 1), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y]).T
        		
                # output array of coefficients
                psfnintp = empty((gdat.numbsidekern, gdat.numbsidekern, kernmatrdesi.shape[1]))

        		# solve p = kernmatrdesi psfnintp for psfnintp
                for iy in gdat.indxsidekern:
                    for ix in gdat.indxsidekern:
        		        p = psf[iy*factkernusam:(iy+1)*factkernusam+1, ix*factkernusam:(ix+1)*factkernusam+1].flatten()
        		        psfnintp[iy, ix, :] = dot(linalg.inv(dot(kernmatrdesi.T, kernmatrdesi)), dot(kernmatrdesi.T, p))
    
        sbrt = dict()
            
        # prepare the evaluation arrays
        if numbtrap > 0:
            if gdat.verbtype > 1:
                print 'Constructing evaluation elements...'
            
            # fill the evaluation dictionary
            initchro(gdat, gdatmodi, strgstat, 'evalelem')
            if strgstat == 'next':
                indxpopleval = gdatmodi.indxpoplmodi
                indxgrideval = array([indxgridpopl[gdatmodi.indxpoplmodi[0]]])
                indxregieval = gdatmodi.indxregimodi
                numbelemeval = gdatmodi.numbelemeval
                indxelemeval = [[[]]]
                print 'gdatmodi.numbelemeval[0][0]'
                print gdatmodi.numbelemeval[0][0]
                indxelemeval[0][0] = arange(gdatmodi.numbelemeval[0][0])
                if not gdatmodi.propfixp:
                    dicteval = gdatmodi.dicteval
                    if gdat.verbtype > 1:
                        for ll, l in enumerate(indxpopleval):
                            for d in indxregipopl[l]:
                                for strgfeat in liststrgfeateval[l]:
                                    print strgfeat
                                    print dicteval[ll][d][strgfeat]
            else:
                indxgrideval = gdat.indxgrid
                indxregieval = gdat.indxregi
                indxpopleval = indxpopl
                indxelemeval = [[[] for d in indxregipopl[l]] for l in indxpopl]
                numbelemeval = numbelem
                # common dictionary
                dicteval = [[{} for d in indxregipopl[l]] for l in indxpopl]
    
                for l in indxpopl:
                    for d in indxregipopl[l]:
                        for strgfeat in liststrgfeateval[l]:
                            if strgfeat == 'spec':
                                indxelemeval[l][d] = arange(numbelem[l][d])
                                if boolelemlghtspat[l]:
                                    dicteval[l][d]['spec'] = retr_spec(gdat, dictelem[l][d]['flux'], sind=dictelem[l][d]['sind'], curv=dictelem[l][d]['curv'], \
                                              expc=dictelem[l][d]['expc'], sind0001=dictelem[l][d]['sind0001'], sind0002=dictelem[l][d]['sind0002'], spectype=spectype[l])
                                if elemtype[l] == 'lghtline':
                                    dicteval[l][d]['spec'] = retr_spec(gdat, dictelem[l][d]['flux'], elin=dictelem[l][d]['elin'], edis=gdat.edis, spectype=spectype[l])
                            else:
                                dicteval[l][d][strgfeat] = dictelem[l][d][strgfeat]
                if gdat.verbtype > 1:
                    for l in indxpopl:
                        for d in indxregipopl[l]:
                            print 'ld'
                            print l, d
                            for strgfeat in liststrgfeateval[l]:
                                print strgfeat
                                print dicteval[l][d][strgfeat]
            stopchro(gdat, gdatmodi, strgstat, 'evalelem')
        else:
            if strgstat == 'next':
                if not gdatmodi.proplenp:
                    indxgrideval = [0]
                if gdatmodi.prophost:
                    indxgrideval = [0, 1]
                if gdatmodi.propsour:
                    indxgrideval = [0, 2]
            else:
                indxgrideval = gdat.indxgrid
                indxregieval = gdat.indxregi
            #numbelemeval = fittnumbelemzero
        
        if strgstat == 'this' or strgstat == 'next' and not gdatmodi.evalllikpert:
            indxpixleval = gdat.indxpixleval
            indxcubeeval = gdat.indxcubeeval
        else:
            indxcubeeval = [[[] for d in indxregieval] for y in indxgrideval]
            # list of pixels for each grid and region
            indxpixleval = [[[] for d in indxregieval] for y in indxgrideval]
        
        # determine the indices of the pixels over which element kernels will be evaluated
        if numbtrap > 0:
            # list of pixels for each grid, population, region, and element
            listindxpixleval = [[[] for d in indxregipopl[l]] for l in indxpopleval]
            print 'listindxpixleval'
            print listindxpixleval
            print 'indxpopleval'
            print indxpopleval
            print 'indxpoplgrid'
            print indxpoplgrid
            print 
            print 'hey'

            # list of pixels for each grid, region, and population -- to be concatenated over populations to obtain indxpixleval
            #indxpixlevalpopl = [[[[] for l in indxpoplgrid[y]] for d in indxregieval] for yy, y in enumerate(indxgrideval)]
            indxpixlevalpopl = [[[] for l in indxpopleval] for d in indxregieval]
            for ll, l in enumerate(indxpopleval):
                for dd, d in enumerate(indxregieval):
                    #for yy, y in enumerate([indxgridpopl[l]]):
                        
                        #print 'yy, ll, l, dd, d'
                        #print yy, ll, l, dd, d
                        #print 'indxpopleval'
                        #print indxpopleval
                        #print 'indxregieval'
                        #print indxregieval
                        #print 'indxgridpopl'
                        #print indxgridpopl
                        #print
                        #print 'numbelemeval'
                        #print numbelemeval
                        #print 'indxregipopl'
                        #print indxregipopl
                        #print
                        #print 'listindxpixleval'
                        #print listindxpixleval
                        #print 'listindxpixleval[yy]'
                        #print listindxpixleval[yy]
                        #print
                        #print 'listindxpixleval[yy][ll]'
                        #print listindxpixleval[yy][ll]
                        #print 'listindxpixleval[yy][ll][dd]'
                        #print listindxpixleval[yy][ll][dd]
                        #print 'indxpixlevalpopl[yy][dd][ll]'
                        #print indxpixlevalpopl[yy][dd][ll]
                        #print
                        #print
                        #print
                        #print
                        
                    if numbelemeval[ll][dd] > 0:
                        if elemspatevaltype[l] == 'locl' or elemspatevaltype[l] == 'loclhash':
                            listindxpixleval[ll][dd], indxpixlevalpopl[dd][ll] = retr_indxpixlevalconc(gdat, strgmodl, dicteval, l, ll, dd)
                if not (strgstat == 'this' or strgstat == 'next' and not gdatmodi.evalllikpert): 
                    print 
                    print
                    print 'indxpixlevalpopl[dd]'
                    print indxpixlevalpopl[dd]
                    print 'dd'
                    print dd
                    print 'indxpixleval'
                    print indxpixleval
                    print 'l'
                    print l
                    print 'indxgridpopl'
                    print indxgridpopl
                    for yy, y in enumerate([indxgridpopl[l]]):
                        indxpixleval[yy][dd] = unique(concatenate(indxpixlevalpopl[dd])).astype(int)
                        indxcubeeval[yy][dd] = meshgrid(indxenereval, indxpixleval[yy][dd], indxevtteval, indexing='ij')
        
                    #    if numbelemeval[ll][dd] > 0:
                    #        if elemspatevaltype[l] == 'locl' or elemspatevaltype[l] == 'loclhash':
                    #            listindxpixleval[yy][ll][dd], indxpixlevalpopl[yy][dd][ll] = retr_indxpixlevalconc(gdat, strgmodl, dicteval, l, ll, dd)
                    #if not (strgstat == 'this' or strgstat == 'next' and not gdatmodi.evalllikpert): 
                    #    indxpixleval[yy][dd] = unique(concatenate(indxpixlevalpopl[yy][dd])).astype(int)
                    #    indxcubeeval[yy][dd] = meshgrid(indxenereval, indxpixleval[yy][dd], indxevtteval, indexing='ij')
        
        # load indxcubeeval to the global object for a possible state update
        if strgstat == 'next':
            gdatmodi.indxregieval = indxregieval
            gdatmodi.indxcubeeval = indxcubeeval
        
        if numbtrap > 0:
            if calcsbrtpntsfrdg:
                sbrtpntsfrgd = [[] for d in enumerate(indxregieval)]
            if calcsbrtextsbgrd:
                sbrtextsbgrd = [[] for d in enumerate(indxregieval)]
            if calcdeflsubh:
                deflelem = [[] for d in indxregieval]
            # retrieve or initialize state variable
            for dd, d in enumerate(indxregieval):
                if calcsbrtpnts:
                    if strgstat == 'next' and not gdatmodi.proppsfp:
                        sbrtpnts[dd] = empty_like(gdat.expo[d])
                        if gdatmodi.propfixp and not gdatmodi.proppsfp:
                            sbrtpnts[d][indxcubeeval[0][dd]] = getattr(gdatobjt, strgpfixthis + 'sbrtpnts')[d][indxcubeeval[0][dd]]
                        if gdatmodi.propelem:
                            sbrtpnts[d][indxcubeeval[0][dd]] = copy(getattr(gdatobjt, strgpfixthis + 'sbrtpnts')[d][indxcubeeval[0][dd]])
                    else:
                        sbrtpnts[d] = zeros_like(gdat.expo[d])
                if calcsbrtextdbgrd:    
                    if strgstat == 'next':
                        sbrtextsbgrd[dd] = empty_like(gdat.expo[d])
                        if gdatmodi.propfixp:
                            sbrtextsbgrd[d][indxcubeeval[0][dd]] = getattr(gdatobjt, strgpfixthis + 'sbrtpnts')[d][indxcubeeval[0][dd]]
                        if gdatmodi.prop:
                            sbrtextsbgrd[d][indxcubeeval[0][dd]] = copy(getattr(gdatobjt, strgpfixthis + 'sbrtpnts')[d][indxcubeeval[0][dd]])
                    else:
                        sbrtextsbgrd[d] = zeros_like(gdat.expo[d])
                if calcdeflelem:
                    if strgstat == 'next':
                        if gdatmodi.propfixp:
                            deflelem[dd] = getattr(gdatobjt, strgpfixthis + 'deflelem')[d]
                        elif gdatmodi.propelem:
                            deflelem[dd] = copy(getattr(gdatobjt, strgpfixthis + 'deflelem')[d])
                    else:
                        deflelem[dd] = zeros((gdat.numbpixl, 2))
            
            if gdat.verbtype > 1:
                print 'elemspatevaltype'
                print elemspatevaltype
                print
        
            # element kernel evaluation
            initchro(gdat, gdatmodi, strgstat, 'kernelem')
            if calcsbrtpnts:
                if strgstat != 'next' or gdatmodi.propelem:
                    for ll, l in enumerate(indxpopleval):
                        if elemtype[l] != 'lens':
                            for dd, d in enumerate(indxregipopl[l]):
                                for k in range(numbelemeval[ll][d]):
                                    if elemtype[l] == 'lght':
                                        varbevalextd = dicteval[ll][d]['spec'][:, k]
                                    if elemtype[l] == 'clus':
                                        varbevalextd = dicteval[ll][d]['nobj'][None, k]
                                    if elemtype[l] == 'lghtline':
                                        varbevalextd = dicteval[ll][dd]['spec'][:, k]
                                    if gdat.verbtype > 1:
                                        print 'varbevalextd'
                                        print varbevalextd
                                        print
                                    if elemtype[l] == 'lght' or elemtype[l] == 'clus':
                                        sbrtpnts[d][:, listindxpixleval[ll][dd][k], :] += retr_sbrtpnts(gdat, dicteval[ll][dd]['lgal'][k], \
                                                                                       dicteval[ll][dd]['bgal'][k], varbevalextd, psfnintp, oaxitype, listindxpixleval[ll][dd][k])
                                    if elemtype[l] == 'lghtline':
                                        sbrtpnts[d][:, 0, 0] += dicteval[ll][dd]['spec'][:, k]
                                        
                    if gdat.diagmode:   
                        cntptempchec = retr_cntp(gdat, sbrtpnts, indxregieval, gdat.indxcubeeval)
                        for d in indxregieval:
                            if amin(cntptempchec[d]) < -0.01:
                                print 'Diagnostic check failed.'
                                print 'cntptempchec'
                                summgene(cntptempchec[d])
                                raise Exception('Element surface brightness is not positive-definite.')
                
                sbrt['pnts'] = [[] for d in indxregieval]
                print 'indxregieval'
                print indxregieval
                for dd, d in enumerate(indxregieval):
                    sbrt['pnts'][dd] = sbrtpnts[d][indxcubeeval[0][dd]]
                
                # when the only background template is the data-PS residual, correct the PS template for numerical noise
                if backtype[0] == 'data':
                    for dd, d in enumerate(indxregieval):
                        sbrt['pnts'][dd][where(sbrt['pnts'][dd] <= 1e-100)] = 1e-100
                setattr(gdatobjt, strgpfix + 'sbrtpnts', sbrt['pnts'])

                if gdat.diagmode:
                    cntppntschec = retr_cntp(gdat, sbrt['pnts'], indxregieval, indxcubeeval)
                    for d in gdat.indxregi:
                        if amin(cntppntschec[d]) < -0.1:
                            print 'Diagnostic check failed.'
                            print 'cntppntschec[d]'
                            summgene(cntppntschec[d])
                            raise Exception('Point source spectral surface brightness is not positive-definite.')
            
            if calcdeflelem:
                if strgstat != 'next' or gdatmodi.propelem:
                    for ll, l in enumerate(indxpopleval):
                        if elemtype[l] == 'lens':
                            for dd, d in enumerate(indxregieval):
                                for kk, k in enumerate(indxelemeval[ll][dd]):
                                    if gdat.variasca:
                                        asca = dicteval[ll][dd]['asca'][k]
                                    else:
                                        asca = gdat.ascaglob
                                    if gdat.variacut:
                                        acut = dicteval[ll][dd]['acut'][k]
                                    else:
                                        acut = gdat.acutglob
                                    
                                    deflelem[dd][listindxpixleval[ll][dd][kk], :] += retr_defl_jitt(gdat.indxpixl, gdat.lgalgrid, gdat.bgalgrid, \
                                                                 dicteval[ll][dd]['lgal'][kk], dicteval[ll][dd]['bgal'][kk], dicteval[ll][dd]['defs'][kk], 0., 0., \
                                                                                                                 asca=asca, acut=acut, indxpixleval=listindxpixleval[ll][dd][kk])
                    setattr(gdatobjt, strgpfix + 'deflelem', deflelem)
                    
                    if gdat.diagmode:
                        for dd, d in enumerate(indxregieval):
                            if not isfinite(deflelem[dd]).all():
                                raise Exception('Element deflection is not finite.')

                for dd, d in enumerate(indxregieval):
                    defl[dd] += deflelem[dd]

            stopchro(gdat, gdatmodi, strgstat, 'kernelem')
        
        sbrtdiff = dict()
        if lensmodltype != 'none':
            
            # lensed surface brightness
            initchro(gdat, gdatmodi, strgstat, 'sbrtlens')
            
            sbrtdiff['lens'] = [[] for d in indxregieval]
            for dd, d in enumerate(indxregieval):
                if gdat.numbener > 1:
                    specsour = retr_spec(gdat, array([fluxsour[d]]), sind=array([sindsour[d]]))
                else:
                    specsour = array([fluxsour[d]])
                sbrtdiff['galxbgrd'][dd] = retr_sbrtsers(gdat, gdat.lgalgrid[indxpixleval[0][dd]], gdat.bgalgrid[indxpixleval[0][dd]], \
                                                                                    lgalsour[d], bgalsour[d], specsour, sizesour[d], ellpsour[d], anglsour[d])
                
                sbrtdiff['bgrd'][dd] = sbrtdiff['galxbgrd'][dd] + sbrtdiff['extsbgrd'][dd]
                sbrtbgrdobjt = sp.interpolate.RectBivariateSpline(gdat.bgalgrid, gdat.lgalgrid, sbrtdiff['bgrd'][dd])
                sbrtdiff['lens'][dd] = sbrtbgrdobjt(gdat.bgalgrid[indxpixleval[0][dd]] - defl[d][indxpixleval[0][dd], 1], \
                                                    gdat.lgalgrid[indxpixleval[0][dd]] - defl[d][indxpixleval[0][dd], 0])

                #sbrtdiff['lens'][dd] = retr_sbrtsers(gdat, gdat.lgalgrid[indxpixleval[0][dd]] - defl[d][indxpixleval[0][dd], 0], \
                #                                           gdat.bgalgrid[indxpixleval[0][dd]] - defl[d][indxpixleval[0][dd], 1], \
                #                                                                    lgalsour[d], bgalsour[d], specsour, sizesour[d], ellpsour[d], anglsour[d])
           
            if gdat.diagmode:
                if not isfinite(sbrtdiff['lens']).all():
                    raise Exception('Lensed emission is not finite.')

            stopchro(gdat, gdatmodi, strgstat, 'sbrtlens')
            setattr(gdatobjt, strgpfix + 'defl', defl)
            
        ### background surface brightness
        numbback = getattr(gdat, strgmodl + 'numbback')
        indxback = getattr(gdat, strgmodl + 'indxback')

        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        indxbacpback = getattr(gdat, strgmodl + 'indxbacpback')
        # temp
        #smthback = getattr(gdat, strgmodl + 'smthback')
        specback = getattr(gdat, strgmodl + 'specback')
        unifback = getattr(gdat, strgmodl + 'unifback')
        sbrtback = [[] for d in range(indxregieval.size)]
        # temp
        #for d in range(indxregieval.size):
        #    sbrtback[d] = empty((numbback, gdat.numbener, indxpixleval[yy][d].size, gdat.numbevtt))
        
        # evaluate host galaxy surface brightness
        if hostemistype != 'none':
            initchro(gdat, gdatmodi, strgstat, 'sbrthost')
            sbrtdiff['host'] = [[] for d in indxregieval]
            for dd, d in enumerate(indxregieval):
                if strgstat == 'next' and not gdatmodi.prophost:
                    if gdat.verbtype > 1:
                        print 'Retrieving the host galaxy surface brightness...'
                    sbrtdiff['host'][dd] = getattr(gdatobjt, strgpfixthis + 'sbrthost')[d][indxcubeeval[0][dd]]
                else:
                    if gdat.verbtype > 1:
                        print 'Host galaxy surface brightness evaluation...'
                    for dd, d in enumerate(indxregieval):
                        if gdat.numbener > 1:
                            spechost = retr_spec(gdat, array([fluxhost[d]]), sind=array([sindhost[d]]))
                        else:
                            spechost = array([fluxhost[d]])
                        
                        if gdat.verbtype > 1:
                            print 'dd, d'
                            print dd, d
                            print 'lgalhost[d]'
                            print lgalhost[d] * gdat.anglfact
                            print 'bgalhost[d]'
                            print bgalhost[d] * gdat.anglfact
                            print 'spechost'
                            print spechost
                            print 'sizehost[d]'
                            print sizehost[d]
                            print 'ellphost[d]'
                            print ellphost[d]
                            print 'anglhost[d]'
                            print anglhost[d]
                            print 'serihost[d]'
                            print serihost[d]
                        sbrtdiff['host'][dd] = retr_sbrtsers(gdat, gdat.lgalgrid, gdat.bgalgrid, lgalhost[d], \
                                                                                    bgalhost[d], spechost, sizehost[d], ellphost[d], anglhost[d], serihost[d])
                    
                    setattr(gdatobjt, strgpfix + 'sbrthost', sbrtdiff['host'])
                if gdat.verbtype > 1:
                    for dd, d in enumerate(indxregieval):
                        print 'sbrtdiff[host][dd]'
                        summgene(sbrtdiff['host'][dd])
            stopchro(gdat, gdatmodi, strgstat, 'sbrthost')
        
        # construct the model diffuse surface brightness
        # convolve the model surface brightness with the PSF
        if gdat.verbtype > 1:
            print 'PSF convolution of diffuse components...'

        initchro(gdat, gdatmodi, strgstat, 'sbrtdiffconv')
        sbrtdiffconv = dict()
        for k, name in enumerate(listnamediff):
            
            if gdat.verbtype > 1:
                print name
           
            booltemp = False
            if strgstat == 'next':
                if name.startswith('back'):
                    if not gdatmodi.proppsfp or unifback[int(name[4:])]:
                        booltemp = True
                elif name == 'host':
                    if not gdatmodi.prophost and not gdatmodi.proppsfp:
                        booltemp = True
            
            if booltemp:
                sbrtdiffconv[name] = [[] for d in indxregieval]
                for dd, d in enumerate(indxregieval):
                    sbrtdiffconv[name][dd] = copy(getattr(gdatobjt, strgpfixthis + 'sbrt' + name + 'conv')[d][indxcubeeval[0][dd]])
                if gdat.verbtype > 1:
                    print 'Retrieving the already-convolved diffuse component from the state vector...'
            else:
                
                if name.startswith('back'):
                    indxbacktemp = int(name[4:8])
                    
                    sbrtdiff[name] = [[] for d in indxregieval]
                    for dd, d in enumerate(indxregieval):
                        if gdat.pixltype == 'heal' and (psfnevaltype == 'full' or psfnevaltype == 'conv') and not unifback[indxbacktemp]:
                            sbrtdiff[name][dd] = getattr(gdat, strgmodl + 'sbrtbackhealfull')[indxbacktemp][d][indxcubeeval[0][dd]]
                        else:
                            sbrtdiff[name][dd] = sbrtbacknorm[indxbacktemp][dd][indxcubeeval[0][dd]]
                
                if gdat.diagmode:
                    if not isfinite(sbrtbacknorm).all():
                        raise Exception('')
                    if not isfinite(sbrtdiff[name]).all():
                        print 'name'
                        print name
                        raise Exception('Diffuse emission is not finite.')
                
                if not (strgstat == 'next' and gdatmodi.propelem and indxgridpopl[gdatmodi.indxpoplmodi[0]] == 0 or not convdiff[k]):
                    sbrtdiffconv[name] = [[] for d in indxregieval]
                    for dd, d in enumerate(indxregieval):
                        sbrtdiffconv[name][dd] = empty_like(indxcubeeval[0][dd][0])
                   
                    # temp -- change this to modified energy and event class bins
                    if gdat.pixltype == 'heal':
                        psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
                        fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
                    for dd, d in enumerate(indxregieval):
                        for ii, i in enumerate(indxenereval):
                            for mm, m in enumerate(indxevtteval):
                                if gdat.pixltype == 'cart':
                                    sbrtdiffconv[name][dd][ii, :, mm] = convolve_fft(sbrtdiff[name][dd][ii, :, mm].reshape((gdat.numbsidecart, gdat.numbsidecart)), \
                                                                                                                                                    psfnconv[m][i]).flatten()
                                    # temp -- potentially slow
                                    indx = where(sbrtdiffconv[name][dd][ii, :, mm] < 1e-50)
                                    sbrtdiffconv[name][dd][ii, indx, mm] = 1e-50
                                if gdat.pixltype == 'heal':
                                    sbrtdiffconv[name][dd][ii, :, mm] = hp.smoothing(sbrtdiff[name][dd][ii, :, mm], fwhm=fwhm[i, m])[gdat.indxpixlrofi]
                                    sbrtdiffconv[name][dd][ii, :, mm][where(sbrtdiffconv[name][dd][ii, :, mm] <= 1e-50)] = 1e-50
                    if gdat.verbtype > 1:
                        print 'Convolving...'
                else:
                    sbrtdiffconv[name] = sbrtdiff[name]
                    if gdat.verbtype > 1:
                        print 'Falling back on the unconvolved components...'
                
                setattr(gdatobjt, strgpfix + 'sbrt' + name + 'conv', sbrtdiffconv[name])
                if gdat.verbtype > 1:
                    for dd, d in enumerate(indxregieval):
                        print 'dd'
                        print dd
                        print 'sbrtdiff[name][dd]'
                        summgene(sbrtdiff[name][dd])
                        print
            if gdat.verbtype > 1:
                for dd, d in enumerate(indxregieval):
                    print 'dd'
                    print dd
                    print 'sbrtdiffconv[name][dd]'
                    summgene(sbrtdiffconv[name][dd])
                    print
        
        stopchro(gdat, gdatmodi, strgstat, 'sbrtdiffconv')
       
        ## model emission
        initchro(gdat, gdatmodi, strgstat, 'sbrtmodl')
        if gdat.verbtype > 1:
            print 'Summing up the model emission...'
        sbrt['modl'] = [[] for d in indxregieval]
        for dd, d in enumerate(indxregieval):
            sbrt['modl'][dd] = zeros_like(indxcubeeval[0][dd][0], dtype=float)
            
        for name in listnamediff:
            sbrt[name] = [[] for d in indxregieval]
            for dd, d in enumerate(indxregieval):
                if name.startswith('back'):
                    indxbacktemp = int(name[4:8])
                    if False:
                        print
                        print
                        print 'indxbacpback'
                        print indxbacpback
                        print 'indxbacktemp'
                        print indxbacktemp
                        print 'indxbacpback[indxbacktemp]'
                        summgene(indxbacpback[indxbacktemp])
                        print 'sbrtdiffconv[name]'
                        summgene(sbrtdiffconv[name])
                        print 'specback[indxbacktemp]'
                        print specback[indxbacktemp]
                        print 'indxbacpback[indxbacktemp]'
                        print indxbacpback[indxbacktemp]
                        print 'indxbacpback[indxbacktemp][indxregieval]'
                        print indxbacpback[indxbacktemp][indxregieval]
                        print 'bacp[indxbacpback[indxbacktemp][indxregieval]]'
                        print bacp[indxbacpback[indxbacktemp][indxregieval]]
                    if specback[indxbacktemp]:
                        sbrt[name][dd] = sbrtdiffconv[name][dd] * bacp[indxbacpback[indxbacktemp][d]]
                    else:
                        if gdat.verbtype > 1:
                            print 'bacp'
                            print bacp
                        sbrt[name][dd] = sbrtdiffconv[name][dd] * bacp[indxbacpback[indxbacktemp]].flatten()[:, None, None]
                else:
                    sbrt[name] = sbrtdiffconv[name]

                sbrt['modl'][dd] += sbrt[name][dd]
                if gdat.verbtype > 1:
                    print name
                    print 'sbrtdiffconv[name][dd]'
                    summgene(sbrtdiffconv[name][dd])
                    print 'sbrt[name][dd]'
                    summgene(sbrt[name][dd])
                    print
        
        if gdat.verbtype > 1:
            for d in gdat.indxregi:
                print 'sbrt[modl][d]'
                summgene(sbrt['modl'][d])
                print

        ## add point source emission
        if numbtrap > 0 and ('lght' in elemtype or 'clus' in elemtype or 'lghtline' in elemtype) and backtype[0] != 'data':
            if gdat.verbtype > 1:
                for d in gdat.indxregi:
                    print 'sbrt[pnts][d]'
                    summgene(sbrt['pnts'][d])
            sbrt['modl'] += sbrt['pnts']
        stopchro(gdat, gdatmodi, strgstat, 'sbrtmodl')
        
        if gdat.verbtype > 1:
            print 'After adding light emitting elements (pnts and line)...'
            for d in gdat.indxregi:
                print 'sbrt[modl][d]'
                summgene(sbrt['modl'][d])
            print

        # temp
        if isinstance(backtype[0], str) and backtype[0].startswith('bfun'):
            for dd, d in enumerate(indxregieval):
                sbrt['modl'][dd][where(sbrt['modl'][dd] < 1e-50)] = 1e-50
        
        ### count map
        initchro(gdat, gdatmodi, strgstat, 'expo')
        cntp = dict()
        cntp['modl'] = retr_cntp(gdat, sbrt['modl'], indxregieval, indxcubeeval)
        setattr(gdatobjt, strgpfix + 'cntpmodl', cntp['modl'])
        stopchro(gdat, gdatmodi, strgstat, 'expo')

        if gdat.verbtype > 1:
            for dd, d in enumerate(indxregieval):
                print 'sbrt[modl][dd]'
                summgene(sbrt['modl'][dd])

        # mock data specific
        if strgmodl == 'true' and strgstat == 'this':
            
            if raww:
                strgvarb = 'truecntpmodlraww'
            else:
                strgvarb = 'cntpdata'
            
            # generate count data
            cntstemp = zeros((gdat.numbregi, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
            for d in gdat.indxregi:
                for i in gdat.indxener:
                    for j in gdat.indxpixl:
                        for m in gdat.indxevtt:
                            cntstemp[d, i, j, m] = poisson(cntp['modl'][d][i, j, m])
            setattr(gdat, strgvarb, cntstemp)
        
            if not gdat.killexpo and amax(cntstemp) == 0 and not raww:
                print 'gdat.deltener'
                summgene(gdat.deltener)
                print 'gdat.expo'
                summgene(gdat.expo)
                print 'gdat.cntpdata'
                summgene(gdat.cntpdata)
                print 'cntp[modl]'
                summgene(cntp['modl'])
                print 'gdat.apix'
                print gdat.apix
                print gdat.apix * gdat.anglfact**2
                raise Exception('Data is zero.')
            
            if raww:
                return
            else:
                retr_datatick(gdat)
       
            proc_cntpdata(gdat)
        
        # diagnostics
        if gdat.diagmode:
            for varb in [sbrt['modl'], cntp['modl']]:
                frac = amin(varb) / mean(varb)
                if frac < -1e-3:
                    print 'frac'
                    print frac
                    summgene(varb)
                    raise Exception('Total model spectral surface brightness is not positive-definite.')

        ### log-likelihood
        initchro(gdat, gdatmodi, strgstat, 'llikcalc')
        #llik = retr_llik_mult(gdat, cntp['modl'])
        
        llik = retr_llik_depr(gdat, cntp['modl'], indxregieval, indxcubeeval)
        
        if gdat.verbtype > 1:
            print 'indxregieval'
            print indxregieval
            print 'indxenereval'
            print indxenereval
            print 'indxevtteval'
            print indxevtteval

        if gdat.verbtype > 1:
            for yy, y in enumerate(indxgrideval):
                for dd, d in enumerate(indxregieval):
                    print 'cntp[modl][dd]'
                    summgene(cntp['modl'][dd])
                    print 'indxpixleval[yy][dd]'
                    summgene(indxpixleval[yy][dd])

        stopchro(gdat, gdatmodi, strgstat, 'llikcalc')
        
        if gdat.diagmode:
            if not isfinite(llik).all():
                raise Exception('Likelihood is not finite.')
    
        if strgstat == 'next':
            thislliktotl = getattr(gdatobjt, strgpfixthis + 'lliktotl')
            thisllik = getattr(gdatobjt, strgpfixthis + 'llik')
            deltlliktotl = 0
            for dd, d in enumerate(indxregieval):
                deltlliktotl += sum(llik[dd] - thisllik[d][indxcubeeval[0][dd]])
            setattr(gdatobjt, strgpfix + 'deltlliktotl', deltlliktotl)
            lliktotl = deltlliktotl + thislliktotl
        else:
            lliktotl = sum(llik)
        setattr(gdatobjt, strgpfix + 'llik', llik)

        if gdat.verbtype > 1:
            if strgstat == 'next':
                print 'thislliktotl'
                print thislliktotl
                print 'deltlliktotl'
                print deltlliktotl
            for dd, d in enumerate(indxregieval):
                print 'llik[dd]'
                summgene(llik[dd])
            print 'lliktotl'
            print lliktotl
            print
        setattr(gdatobjt, strgpfix + 'lliktotl', lliktotl) 
    else:
        if gdat.verbtype > 1:
            print 'Not evaluating the likelihood...'
        gdatmodi.nextdeltlliktotl = 0.
    
    stopchro(gdat, gdatmodi, strgstat, 'llik')

    # log-prior
    initchro(gdat, gdatmodi, strgstat, 'lpri')
    if gdat.verbtype > 1:
        print 'Evaluating the prior...'
        
    lpri = zeros(gdat.numblpri)
    if numbtrap > 0:
        
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
    
        meanelem = sampvarb[indxfixpmeanelem]
        
        if gdat.penalpridiff:
            
            if strgstat == 'next' and gdatmodi.prophypr:
                lpri[1] = gdatmodi.thislpridiff
            else:
                sbrtdatapnts = gdat.sbrtdata - sbrt['pnts']
                if gdat.pixltype == 'heal':
                    raise Exception('')
                if gdat.pixltype == 'cart':
                    psecodimdatapnts = empty((gdat.numbregi, gdat.numbener, gdat.numbsidecart / 2, gdat.numbevtt))
                    psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
                    fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
                    sigm = fwhm / 2.355
                    psecodimdatapntsprio = exp(-2. * gdat.meanmpolodim[None, :, None] / (0.1 / sigm[:, None, :]))
                    lpridiff = 0.
                    for d in gdat.indxregi:
                        for i in gdat.indxener:
                            for m in gdat.indxevtt:
                                psecdatapnts = retr_psec(gdat, sbrtdatapnts[d, i, :, m])
                                psecodimdatapnts[d, i, :, m] = retr_psecodim(gdat, psecdatapnts)
                                psecodimdatapnts[d, i, :, m] /= psecodimdatapnts[d, i, 0, m]
                                lpridiff += -0.5 * sum((psecodimdatapnts[i, :, m] - psecodimdatapntsprio[d, i, :, m])**2)
                                setattr(gdatobjt, strgpfix + 'psecodimdatapntsreg%dene%devt%d' % (d, i, m), psecodimdatapnts[d, i, :, m])
                                setattr(gdatobjt, strgpfix + 'psecodimdatapntsprioreg%dene%devt%d'% (d, i, m), psecodimdatapntsprio[d, i, :, m])
                lpri[1] = lpridiff 
                if strgstat == 'this' or strgstat == 'next':
                    setattr(gdatobjt, strgpfix + 'lpridiff', lpridiff)
        
        for l in indxpopl:
            for d in indxregipopl[l]:

                lpri[0] -= 0.5 * gdat.priofactdoff * numbcomp[l] * numbelem[l][d]
                lpri[2] += retr_lprbpois(numbelem[l][d], meanelem[l])
                
                for k, (strgfeat, strgpdfn) in enumerate(zip(liststrgfeatprio[l], liststrgpdfnprio[l])):
                    
                    indxlpri = 3 + l * gdat.numbregi * maxmnumbcomp + d * maxmnumbcomp + k
                    if False and strgpdfn == 'tmpl':

                        if strgpdfn.endswith('cons'):
                            pdfnspatpriotemp = getattr(gdat, strgmodl + 'pdfnspatpriotemp')
                            spatdistcons = sampvarb[getattr(gdat, strgmodl + 'indxfixpspatdistcons')]
                            lpdfspatprio, lpdfspatprioobjt = retr_spatprio(gdat, pdfnspatpriotemp, spatdistcons)
                            lpdfspatpriointp = lpdfspatprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
                            
                            # temp
                            lpdfspatpriointp = lpdfspatpriointp.T
                            
                            setattr(gdatobjt, strgpfix + 'lpdfspatpriointp', lpdfspatpriointp)
                            setattr(gdatobjt, strgpfix + 'lpdfspatprioobjt', lpdfspatprioobjt)
                
                        else:
                            lpdfspatprioobjt = gdat.fittlpdfspatprioobjt
                    
                    if strgpdfn == 'self':
                        minmfeat = getattr(gdat, 'minm' + strgfeat)
                        maxmfeat = getattr(gdat, 'maxm' + strgfeat)
                        lpri[indxlpri] = numbelem[l][d] * log(1. / (maxmfeat - minmfeat))
                    elif strgpdfn == 'disc':
                        indxfixpbgaldistscal = getattr(gdat, strgmodl + 'indxfixpbgaldistscalpop%d' % l)
                        lpri[indxlpri] = sum(log(pdfn_dexp(dictelem[l][d]['bgal'], gdat.maxmgang, sampvarb[indxfixpbgaldistscal]))) 
                    elif strgpdfn == 'exposcal':
                        gang = retr_gang(dictelem[l][d]['lgal'], dictelem[l][d]['bgal'])
                        indxfixpgangdistscal = getattr(gdat, strgmodl + 'indxfixpgangdistscalpop%d' % l)
                        lpri[indxlpri] = sum(log(pdfn_expo(gang, gdat.maxmgang, sampvarb[indxfixpgangdistscal]))) 
                        lpri[indxlpri] = -numbelem[l][d] * log(2. * pi) 
                    elif strgpdfn == 'tmpl':
                        lpri[indxlpri] = sum(lpdfspatprioobjt(dictelem[l][d]['bgal'], dictelem[l][d]['lgal'], grid=False))
                    
                    # temp
                    #if gdat.mask != None:
                    #    indxbadd = where((dicttemp['lgalconc'] > gdat.mask[0]) & (dicttemp['lgalconc'] < gdat.mask[1]) & \
                    #                                (dicttemp['bgalconc'] > gdat.mask[2]) & (dicttemp['bgalconc'] < gdat.mask[3]))[0]
                    #    lpri[0] -= 1e6 * indxbadd.size

                    # temp -- this should not be here
                    elif strgpdfn == 'powrslop':
                        lpri[indxlpri] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, dictelem[l][d][strgfeat], strgfeat, sampvarb, l)
                    elif strgpdfn == 'gaus':
                        lpri[indxlpri] = retr_lprigausdist(gdat, gdatmodi, strgmodl, dictelem[l][d][strgfeat], strgfeat, sampvarb, l)
            
        if strgstat == 'next':
            if gdatmodi.propcomp:
                lpridist = 0.
                strgcomp = liststrgcomp[gdatmodi.indxpoplmodi[0]][gdatmodi.indxcompmodi]
                strgpdfn = liststrgpdfnprio[gdatmodi.indxpoplmodi[0]][gdatmodi.indxcompmodi]
                if strgpdfn == 'self':
                    lpridist = retr_lpriselfdist(gdat, strgmodl, gdatmodi.thiscomp, strgcomp) - \
                                                    retr_lpriselfdist(gdat, strgmodl, gdatmodi.nextcomp, strgcomp)
                if strgpdfn == 'powrslop':
                    lpridist = retr_lpripowrdist(gdat, gdatmodi, strgmodl, gdatmodi.thiscomp, strgcomp, sampvarb, gdatmodi.indxpoplmodi) - \
                                                    retr_lpripowrdist(gdat, gdatmodi, strgmodl, gdatmodi.nextcomp, strgcomp, sampvarb, gdatmodi.indxpoplmodi)
                if strgpdfn == 'gaus':
                    lpridist = retr_lprigausdist(gdat, gdatmodi, strgmodl, dictelem[l][d][strgfeat], strgcomp, sampvarb, l)
                if strgpdfn == 'igam':
                    lpridist = retr_lpriigamdist(gdat, gdatmodi, strgmodl, dictelem[l][d][strgfeat], strgcomp, sampvarb, l)
                setattr(gdatobjt, 'nextlpridist', lpridist)
                if gdat.verbtype > 1:
                    print 'lpridist'
                    print lpridist

            if gdatmodi.proptran:
                lpau = zeros(gdat.numblpau) 
                
                thissampvarb = getattr(gdatobjt, strgpfixthis + 'sampvarb')
                if gdatmodi.propbrth or gdatmodi.propsplt:
                    sampvarbtemp = getattr(gdatobjt, strgpfix + 'sampvarb')
                if gdatmodi.propdeth or gdatmodi.propmerg:
                    sampvarbtemp = getattr(gdatobjt, strgpfixthis + 'sampvarb')
        
                if gdatmodi.propbrth or gdatmodi.propdeth:
                
                    for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
                        if listscalcomp[gdatmodi.indxpoplmodi[0]][k] == 'self':
                            minm = getattr(gdat, 'minm' + strgcomp)
                            maxm = getattr(gdat, 'maxm' + strgcomp)
                            lpau[k] = -log(maxm - minm)
                        if listscalcomp[gdatmodi.indxpoplmodi[0]][k] == 'powrslop':
                            lpau[k] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                               strgcomp, thissampvarb, gdatmodi.indxpoplmodi[0])
                        if listscalcomp[gdatmodi.indxpoplmodi[0]][k] == 'gaus':
                            lpau[k] = retr_lprigausdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                               strgcomp, thissampvarb, gdatmodi.indxpoplmodi[0])
                if gdatmodi.propsplt or gdatmodi.propmerg:
                
                    for k, strgcomp in enumerate(liststrgcomp[gdatmodi.indxpoplmodi[0]]):
                        if k == 0 or k == 1:
                            lpau[gdat.fittmaxmnumbcomp*l+k] = -log(2. *pi) - log(gdat.radispmr) -0.5 * (gdatmodi.auxipara[k] / gdat.radispmr)
                        elif k == 2:
                            if gdatmodi.auxipara[k] < 0. or gdatmodi.auxipara[k] > 1.:
                                lpau[k] = -inf
                            else:
                                lpau[k] = 0.
                        elif listscalcomp[gdatmodi.indxpoplmodi[0]][k] == 'powrslop':
                            lpau[k] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                               strgcomp, thissampvarb, gdatmodi.indxpoplmodi[0])
                        elif listscalcomp[gdatmodi.indxpoplmodi[0]][k] == 'gaus':
                            lpau[k] = retr_lprigausdist(gdat, gdatmodi, strgmodl, sampvarbtemp[gdatmodi.indxsamptran[0][k]], \
                                                                                               strgcomp, thissampvarb, gdatmodi.indxpoplmodi[0])
                    
                if gdatmodi.propbrth or gdatmodi.propsplt:
                    lpau *= -1.
                
                lpautotl = sum(lpau)
                if gdat.verbtype > 1:
                    print 'lpau'
                    print lpau
                    print 'lpautotl'
                    print lpautotl

                setattr(gdatobjt, 'nextlpau', lpau)
                setattr(gdatobjt, 'nextlpautotl', lpautotl)
     
    lpritotl = sum(lpri)

    if strgstat == 'next':
        thislpritotl = getattr(gdatobjt, strgpfixthis + 'lpritotl')
        if gdat.verbtype > 1:
            print 'lpri'
            print lpri
            print 'lpritotl'
            print lpritotl
            print 'thislpritotl'
            print thislpritotl
        deltlpritotl = lpritotl - thislpritotl
        if gdatmodi.proptran:
            deltlpritotl += gdatmodi.nextlpautotl
        if gdatmodi.propcomp:
            deltlpritotl += gdatmodi.nextlpridist
        if (gdatmodi.propcomp or gdatmodi.proppsfp or gdatmodi.propbacp or gdatmodi.proplenp) and abs(deltlpritotl) > 1e-3:# or abs(thislpritotl - thissampvarb[0] * 47) > 10.:
            print 'gdatmodi.nextlpridist'
            print gdatmodi.nextlpridist
            print 'deltlpritotl'
            print deltlpritotl
            raise Exception('')
        setattr(gdatobjt, strgpfix + 'deltlpritotl', deltlpritotl)

    setattr(gdatobjt, strgpfix + 'lpritotl', lpritotl) 
    setattr(gdatobjt, strgpfix + 'lpri', lpri)
    
    stopchro(gdat, gdatmodi, strgstat, 'lpri')
    
    if strgstat != 'next':
        lpostotl = lpritotl + lliktotl
        setattr(gdatobjt, strgpfix + 'lpostotl', lpostotl) 
    
    if strgstat == 'next':
        setattr(gdatmodi, 'thislpriprop', lpri)
    
    stopchro(gdat, gdatmodi, strgstat, 'proc')
    
    if fast:
        return
    
    ## tertiary variables that are not needed for evolving the chain
    if strgstat != 'next':
        
        dicttert = [dict() for d in gdat.indxregi]
        
        if strgstat == 'this' and numbtrap > 0:
            numbelempopl = zeros(numbpopl)
            for l in indxpopl:
                numbelempopl[l] = sum(numbelem[l])
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
      
        ## set previously calculated tertiary variables
        setattr(gdatobjt, strgpfix + 'bacp', bacp)
        if psfnevaltype != 'none':
            setattr(gdatobjt, strgpfix + 'psfp', psfp)
        
        ## derived variables
        ## residual count map 
        cntp['resi'] = [[] for d in gdat.indxregi]
        for d in gdat.indxregi:
            cntp['resi'][d] = gdat.cntpdata[d] - cntp['modl'][d]
        setattr(gdatobjt, strgpfix + 'cntpresi', cntp['resi'])
        
        namefeatsort = getattr(gdat, strgmodl + 'namefeatsort')
        
        if lensmodltype != 'none':
            masshostbein = array([gdat.massfrombein * beinhost**2])
            setattr(gdatobjt, strgpfix + 'masshostbein', masshostbein)
            ### sort with respect to deflection at scale radius
            if numbtrap > 0:
                for l in indxpopl:
                    for d in indxregipopl[l]:
                        if numbelem[l][d] > 0:
                            indxelemsortampl = argsort(dictelem[l][d][namefeatsort[l]])[::-1]
                            for strgcomp in liststrgcomp[l]:
                                dictelem[l][d][strgcomp + 'sort'] = dictelem[l][d][strgcomp][indxelemsortampl]

            deflsing = [[] for d in indxregieval]
            for d in indxregieval:
                deflsing[d] = zeros((gdat.numbpixl, 2, numbdeflsingplot))
            conv = zeros((gdat.numbregi, gdat.numbpixl))
            convpsec = zeros((gdat.numbregi, (gdat.numbsidecart / 2)**2))
            convpsecodim = zeros((gdat.numbregi, gdat.numbsidecart / 2))
            if numbtrap > 0:
                indxpopllens = elemtype.index('lens')
            for d in gdat.indxregi:
                numbdeflsing = 2
                if numbtrap > 0:
                    if numbelem[indxpopllens][d] > 0:
                        numbdeflsing += min(numbdeflsubhplot, numbelem[indxpopllens][d]) 
                        numbdeflsing += 1
                for k in range(numbdeflsing):
                    if k == 0:
                        deflsing[d][:, :, k] = deflhost[d]
                    elif k == 1:
                        deflsing[d][:, :, k] = deflextr[d]
                    elif k == 2:
                        deflsing[d][:, :, k] = defl[d] - deflextr[d] - deflhost[d]
                    else:
                        if gdat.variasca:
                            asca = dictelem[indxpopllens][d]['ascasort'][k-3]
                        else:
                            asca = gdat.ascaglob
                        if gdat.variacut:
                            acut = dictelem[indxpopllens][d]['acutsort'][k-3]
                        else:
                            acut = gdat.acutglob
                        deflsing[d][indxpixleval[1][d], :, k] = retr_defl(gdat, dictelem[indxpopllens][d]['lgalsort'][k-3], dictelem[indxpopllens][d]['bgalsort'][k-3], \
                                                                                    dictelem[indxpopllens][d]['defssort'][k-3], 0., 0., \
                                                                                                                   asca=asca, acut=acut, indxpixleval=indxpixleval[1][d])

                # convergence
                ## total
                conv[d, :] = retr_conv(gdat, defl[d]) 
                
                ### power spectrum
                #### two dimensional
                convpsec[d, :] = retr_psec(gdat, conv[d, :])
                
                #### one dimensional
                convpsecodim[d, :] = retr_psecodim(gdat, convpsec[d, :]) 
            setattr(gdatobjt, strgpfix + 'conv', conv)
            setattr(gdatobjt, strgpfix + 'convpsec', convpsec)
            setattr(gdatobjt, strgpfix + 'convpsecodim', convpsecodim)
            
            ## subhalos
            convelem = zeros((gdat.numbregi, gdat.numbpixl))
            convpsecelem = zeros((gdat.numbregi, (gdat.numbsidecart / 2)**2))
            convpsecelemodim = zeros((gdat.numbregi, gdat.numbsidecart / 2))
            if numbtrap > 0:
                for d in gdat.indxregi:
                    ### convergence
                    convelem[d, :] = retr_conv(gdat, deflelem[d]) 
                    
                    ###  power spectrum
                    ##### two dimensional
                    convpsecelem[d, :] = retr_psec(gdat, convelem[d, :])
                    ##### one dimensional
                    convpsecelemodim[d, :] = retr_psecodim(gdat, convpsecelem[d, :]) 
            setattr(gdatobjt, strgpfix + 'convpsecelem', convpsecelem)
            setattr(gdatobjt, strgpfix + 'convpsecelemodim', convpsecelemodim)
            setattr(gdatobjt, strgpfix + 'convelem', convelem)
            
            ### magnification
            magn = empty((gdat.numbregi, gdat.numbpixl))
            histdefl = empty((gdat.numbregi, gdat.numbdefl))
            histdeflelem = empty((gdat.numbregi, gdat.numbdefl))
            deflsingmgtd = zeros((gdat.numbregi, gdat.numbpixl, numbdeflsingplot))
            for d in gdat.indxregi:
                magn[d, :] = 1. / retr_invm(gdat, defl[d]) 
                histdefl[d, :] = histogram(defl[d], bins=gdat.binsdefl)[0]
                if numbtrap > 0:
                    histdeflelem[d, :] = histogram(deflelem[d], bins=gdat.binsdeflelem)[0]
                deflsingmgtd[d, :, :] = sqrt(sum(deflsing[d]**2, axis=1))
            setattr(gdatobjt, strgpfix + 'magn', magn)
            if numbtrap > 0:
                setattr(gdatobjt, strgpfix + 'histdeflelem', histdeflelem)
            setattr(gdatobjt, strgpfix + 'deflsing', deflsing)
            setattr(gdatobjt, strgpfix + 'deflsingmgtd', deflsingmgtd)
            setattr(gdatobjt, strgpfix + 'histdefl', histdefl)
        
        liststrgfeatodim = getattr(gdat, strgmodl + 'liststrgfeatodim')
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
        liststrgfeat = getattr(gdat, strgmodl + 'liststrgfeat')
        
        ## element related
        if numbtrap > 0:
            if False and gdat.numbpixl == 1:
                for l in indxpopl:
                    for d in indxregipopl[l]:
                        for k in range(numbelem[l][d]):
                            setattr(gdatobjt, strgpfix + 'speclinepop%dreg%d%04d' % (l, d, k), dictelem[l][d]['spec'][:, k])
            
            if gdat.datatype == 'mock' and strgmodl == 'true' and gdat.numbpixl > 1:
                gdat.refrlgal = [[[] for d in gdat.trueindxregipopl[l]] for l in gdat.trueindxpopl]
                gdat.refrbgal = [[[] for d in gdat.trueindxregipopl[l]] for l in gdat.trueindxpopl]
                for l in gdat.trueindxpopl:
                    for d in gdat.trueindxregipopl[l]:
                        gdat.refrlgal[l][d] = tile(dictelem[l][d]['lgal'], [3] + list(ones(dictelem[l][d]['lgal'].ndim, dtype=int)))
                        gdat.refrbgal[l][d] = tile(dictelem[l][d]['bgal'], [3] + list(ones(dictelem[l][d]['bgal'].ndim, dtype=int)))
    
            # correlate the fitting model elements with the reference elements
            if gdat.refrinfo and not (strgmodl == 'true' and gdat.datatype == 'mock') and gdat.asscrefr:
                indxelemrefrasschits = [[[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                indxelemfittasschits = [[[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        for d in gdat.fittindxregipopl[l]:
                            if gdat.refrnumbelem[q][d] == 0:
                                continue

                            indxelemfittmatr = empty((gdat.refrnumbelem[q][d], numbelem[l][d]), dtype=int)
                            indxelemrefrmatr = empty((gdat.refrnumbelem[q][d], numbelem[l][d]), dtype=int)
                            matrdist = empty((gdat.refrnumbelem[q][d], numbelem[l][d]))
                            for k in range(numbelem[l][d]):
                                # construct a matrix of angular distances between reference and fitting elements
                                if elemtype[l] == 'lghtline':
                                    matrdist[:, k] = gdat.refrelin[q][d][0, :] - dictelem[l][d]['elin'][k]
                                else:
                                    matrdist[:, k] = retr_angldist(gdat, gdat.refrlgal[q][d][0, :], gdat.refrbgal[q][d][0, :], dictelem[l][d]['lgal'][k], dictelem[l][d]['bgal'][k])
                                indxelemrefrmatr[:, k] = arange(gdat.refrnumbelem[q][d])
                                indxelemfittmatr[:, k] = k
                            matrdist = matrdist.flatten()
                            indxelemrefrmatr = indxelemrefrmatr.flatten()
                            indxelemfittmatr = indxelemfittmatr.flatten()

                            # take only angular separations smaller than some threshold
                            indxmatrthrs = where(matrdist < gdat.anglassc)
                            matrdist = matrdist[indxmatrthrs]
                            indxelemrefrmatr = indxelemrefrmatr[indxmatrthrs]
                            indxelemfittmatr = indxelemfittmatr[indxmatrthrs]

                            # sort the remaining associations with respect to distance
                            indxmatrsort = argsort(matrdist)
                            matrdist = matrdist[indxmatrsort]
                            indxelemrefrmatr = indxelemrefrmatr[indxmatrsort]
                            indxelemfittmatr = indxelemfittmatr[indxmatrsort]
                            
                            for c in range(matrdist.size):
                                if indxelemrefrmatr[c] in indxelemrefrasschits[q][l][d] or indxelemfittmatr[c] in indxelemfittasschits[q][l][d]:
                                    continue
                                indxelemrefrasschits[q][l][d].append(indxelemrefrmatr[c])
                                indxelemfittasschits[q][l][d].append(indxelemfittmatr[c])
                            
                            indxelemrefrasschits[q][l][d] = array(indxelemrefrasschits[q][l][d])
                            indxelemfittasschits[q][l][d] = array(indxelemfittasschits[q][l][d])
                setattr(gdatobjt, strgpfix + 'indxelemrefrasschits', indxelemrefrasschits)
                setattr(gdatobjt, strgpfix + 'indxelemfittasschits', indxelemfittasschits)
                
                indxelemrefrasscmiss = [[[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                indxelemfittasscfals = [[[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        for d in gdat.fittindxregipopl[l]:
                            # indices of the reference elements not associated with the fitting model elements
                            if gdat.refrnumbelem[q][d] > 0:
                                indxelemrefrasscmiss[q][l][d] = setdiff1d(arange(gdat.refrnumbelem[q][d]), indxelemrefrasschits[q][l][d])
                            # indices of the fitting model elements not associated with the reference elements
                            if numbelem[l][d] > 0:
                                indxelemfittasscfals[q][l][d] = setdiff1d(arange(numbelem[l][d]), indxelemfittasschits[q][l][d])
                setattr(gdatobjt, strgpfix + 'indxelemrefrasscmiss', indxelemrefrasscmiss)
                setattr(gdatobjt, strgpfix + 'indxelemfittasscfals', indxelemfittasscfals)
                
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        # collect the associated reference element feature for each fitting element 
                        for namefeatrefr in gdat.refrliststrgfeatonly[q][l]:
                            name = namefeatrefr + gdat.listnamerefr[q]
                            refrfeat = getattr(gdat, 'refr' + namefeatrefr)
                            for d in gdat.fittindxregipopl[l]:
                                dictelem[l][d][name] = zeros(numbelem[l][d])
                                if len(refrfeat[q][d]) > 0 and len(indxelemrefrasschits[q][l][d]) > 0:
                                    dictelem[l][d][name][indxelemfittasschits[q][l][d]] = refrfeat[q][d][0, indxelemrefrasschits[q][l][d]]

            # region indices for elements
            indxregielem = [[[] for d in indxregipopl[l]] for l in indxpopl]
            for l in indxpopl:
                for d in indxregipopl[l]:
                    if False:
                        print 'ld'
                        print l, d
                        print 'numbelem'
                        print numbelem
                        print 'full([numbelem[l][d]], d, dtype=int)'
                        print full([numbelem[l][d]], d, dtype=int)
                    indxregielem[l][d] = full([numbelem[l][d]], d, dtype=int)
            
            # luminosity
            for l in indxpopl:
                if 'reds' in liststrgfeat[l] and 'flux' in liststrgfeat[l]:
                    for d in indxregipopl[l]:
                        ldis = (1. + dictelem[l][d]['redswo08'])**2 * gdat.adisobjt(dictelem[l][d]['redswo08'])
                        lumi = 4. * pi * ldis**2 * dictelem[l][d]['flux']
                        dictelem[l][d]['lumi'] = lumi

            ### derived quantities
            for l in indxpopl:
                for d in indxregipopl[l]:
                    if gdat.numbpixl > 1:
                        #### radial and angular coordinates
                        dictelem[l][d]['gang'] = retr_gang(dictelem[l][d]['lgal'], dictelem[l][d]['bgal'])
                        dictelem[l][d]['aang'] = retr_aang(dictelem[l][d]['lgal'], dictelem[l][d]['bgal'])
                   
                    if elemtype[l] == 'lght':
                        #### number of expected counts
                        dictelem[l][d]['cnts'] = retr_cntspnts(gdat, dictelem[l][d]['lgal'], dictelem[l][d]['bgal'], dictelem[l][d]['spec'], indxregielem[l][d])
                
                #### delta log-likelihood
                for l in indxpopl:
                    dictelem[l][d]['deltllik'] = zeros(numbelem[l][d])
                if gdat.calcllik and not (strgmodl == 'true' and gdat.checprio): 
                    if gdat.verbtype > 1:
                        print
                        print 'Calculating log-likelihood differences when removing elements from the model.'
                    gdatmoditemp = tdpy.util.gdatstrt()
                    for l in indxpopl:
                        if gdat.verbtype > 1:
                            print 'l'
                            print l
                        for k in range(numbelem[l][d]):
                            prep_gdatmodi(gdat, gdatmoditemp, gdatobjt, strgstat, strgmodl)
                            prop_stat(gdat, gdatmoditemp, strgmodl, deth=True, thisindxpopl=l, thisindxregi=d, thisindxelem=k)
                            
                            proc_samp(gdat, gdatmoditemp, 'next', strgmodl)
                            if gdat.diagmode:
                                if not isfinite(lliktotl):
                                    print 'lliktotl'
                                    print lliktotl
                                    raise Exception('')
                                if strgstat == 'next' and gdatmodi.propllik and not isfinite(gdatmoditemp.nextlliktotl):
                                    print 'gdatmoditemp.nextlliktotl'
                                    print gdatmoditemp.nextlliktotl
                                    raise Exception('')
                            
                            if strgmodl == 'true':
                                nextlliktotl = gdat.nexttruelliktotl
                            else:
                                nextlliktotl = gdatmoditemp.nextlliktotl
                            dictelem[l][d]['deltllik'][k] = lliktotl - nextlliktotl
                            if gdat.verbtype > 1:
                                print 'k'
                                print k
                                print 'nextlliktotl'
                                print nextlliktotl
                                print 'lliktotl'
                                print lliktotl
                                print
                        if gdat.verbtype > 1:
                            print
                    if gdat.verbtype > 1:
                        print 'deltllik calculation ended.'
                        print
        
        calcerrr = getattr(gdat, strgmodl + 'calcerrr')

        # more derived parameters
        if numbtrap > 0 and ('lght' in elemtype or 'clus' in elemtype):
            if oaxitype:
                setattr(gdatobjt, strgpfix + 'factoaxi', factoaxi)
           
            setattr(gdatobjt, strgpfix + 'psfn', psfn)

            ### PSF FWHM
            fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            setattr(gdatobjt, strgpfix + 'fwhm', fwhm)
            
            if gdat.fittnumbtrap > 0:
                sbrt['pntssubt'] = zeros_like(gdat.expo)
                for l in indxpopl:
                    for d in indxregipopl[l]:
                        dicteval[l][d]['deltllik'] = dictelem[l][d]['deltllik']
                    if calcerrr[l]:
                        numbiter = 2
                    else:
                        numbiter = 1
                    if calcerrr[l]:
                        sbrt['pntsfull'] = zeros_like(gdat.expo)
                    for d in indxregipopl[l]:
                        for a in range(numbiter):
                            for k in range(numbelemeval[l][d]):
                                if elemtype[l] == 'lght':
                                    varbevalextd = dicteval[l][d]['spec'][:, k]
                                if elemtype[l] == 'clus':
                                    varbevalextd = dicteval[l][d]['nobj'][None, k]
                                if a == 0:
                                    # temp
                                    if dicteval[l][d]['deltllik'][k] < 35:
                                        sbrt['pntssubt'][d][:, listindxpixleval[l][d][k], :] += retr_sbrtpnts(gdat, dicteval[l][d]['lgal'][k], dicteval[l][d]['bgal'][k], \
                                                                                                                varbevalextd, psfnintp, oaxitype, listindxpixleval[l][d][k])
                                if a == 1:
                                    sbrt['pntsfull'][d][:, :, :] += retr_sbrtpnts(gdat, dicteval[l][d]['lgal'][k], dicteval[l][d]['bgal'][k], \
                                                                                                                    varbevalextd, psfnintp, oaxitype, gdat.indxpixl)
                
                        setattr(gdatobjt, strgpfix + 'sbrtpntssubtreg%dpop%d' % (d, l), sbrt['pntssubt'][d])
                    if calcerrr[l]:
                        cntppntsfull = retr_cntp(gdat, sbrt['pntsfull'], gdat.indxregi, gdat.indxcubeeval)
                        cntppnts = retr_cntp(gdat, sbrt['pnts'], indxregieval, indxcubeeval)
                        cntperrr = [[] for d in gdat.indxregi]
                        for d in gdat.indxregi:
                            cntperrr[d] = cntppnts[d] - cntppntsfull[d]
                            if amax(fabs(cntperrr[d] / cntppntsfull[d])) > 0.1:
                                raise Exception('Approximation error in calculating the PS surface brightness is above the tolerance level.')
                            setattr(gdatobjt, strgpfix + 'cntperrrreg%dpop%d' % (d, l), cntperrr[d])

            #fluxbrgt, fluxbrgtassc = retr_fluxbrgt(gdat, dictconc['lgalconc'], dictconc['bgalconc'], concatenate(dictconc['flux']))
            #setattr(gdatobjt, strgpfix + 'fluxbrgt', fluxbrgt)
            #setattr(gdatobjt, strgpfix + 'fluxbrgtassc', fluxbrgtassc)
        
        #### count maps
        if gdat.correxpo:
            cntp = dict()
            for nameecom in listnameecomtotl:
                cntp[nameecom] = retr_cntp(gdat, sbrt[nameecom], indxregieval, indxcubeeval)
                setattr(gdatobjt, strgpfix + 'cntp' + nameecom, cntp[nameecom])
        
        ### spatial averages
        sbrtmean = dict()
        for nameecom in listnameecomtotl:
            sbrtmean[nameecom] = retr_spatmean(gdat, sbrt[nameecom], indxregieval)
            setattr(gdatobjt, strgpfix + 'sbrt%smean' % nameecom, sbrtmean[nameecom])

        # find the 1-point function of the count maps of all emission components including the total emission
        for nameecom in listnameecomtotl:
            for d in gdat.indxregi: 
                for m in gdat.indxevtt:
                    if gdat.numbpixl > 1:
                        for i in gdat.indxener: 
                            histcntp = histogram(cntp[nameecom][d][i, :, m], bins=gdat.binscntpdata)[0]
                            setattr(gdatobjt, strgpfix + 'histcntp' + nameecom + 'reg%dene%devt%d' % (d, i, m), histcntp)
                    else:
                        histcntp = histogram(cntp[nameecom][d][:, 0, m], bins=gdat.binscntpdata)[0]
                        setattr(gdatobjt, strgpfix + 'histcntp' + nameecom + 'reg%devt%d' % (d, m), histcntp)

        if lensmodltype != 'none':
            if strgmodl == 'true':
                s2nr = cntp['lens'] / sqrt(cntp['modl'])
                setattr(gdatobjt, strgpfix + 's2nr', s2nr)
            deflmgtd = sqrt(sum(defl**2, axis=2))
            cntplensgrad = empty((gdat.numbregi, gdat.numbener, gdat.numbpixl, gdat.numbevtt, 2))
            for d in gdat.indxregi:
                for i in gdat.indxener:
                    for m in gdat.indxevtt:
                        cntplensgrad[d, i, :, m, :] = retr_gradmaps(gdat, cntp['lens'][d][i, :, m]) * gdat.sizepixl
            cntplensgradmgtd = sqrt(sum(cntplensgrad**2, axis=4))
            cntplensgrad *= gdat.sizepixl
            indx = where(fabs(cntplensgrad) > 1. * gdat.sizepixl)
            cntplensgrad[indx] = sign(cntplensgrad[indx]) * 1. * gdat.sizepixl
            setattr(gdatobjt, strgpfix + 'deflmgtd', deflmgtd)
            setattr(gdatobjt, strgpfix + 'cntplensgrad', cntplensgrad)
            setattr(gdatobjt, strgpfix + 'cntplensgradmgtd', cntplensgradmgtd)

        if numbtrap > 0:
            for l in indxpopl:
                for d in indxregipopl[l]:
                    if elemtype[l] == 'lght':
                        #### spectra
                        dictelem[l][d]['specplot'] = retr_spec(gdat, dictelem[l][d]['flux'], sind=dictelem[l][d]['sind'], \
                                                                         curv=dictelem[l][d]['curv'], expc=dictelem[l][d]['expc'], \
                                                                         sind0001=dictelem[l][d]['sind0001'], sind0002=dictelem[l][d]['sind0002'], spectype=spectype[l], plot=True)
                        
                        if gdat.datatype == 'inpt':
                            if gdat.exprtype == 'ferm':
                                if False:
                                    print 'dictelem[l][d][lgal]'
                                    print dictelem[l][d]['lgal']
                                    print 'dictelem[l][d][bgal]'
                                    print dictelem[l][d]['bgal']
                                    print 'gdat.sbrt0018objt(dictelem[l][d][bgal], dictelem[l][d][lgal])'
                                    print gdat.sbrt0018objt(dictelem[l][d]['bgal'], dictelem[l][d]['lgal'])
                                    print ''
                                # temp
                                try:
                                    dictelem[l][d]['sbrt0018'] = gdat.sbrt0018objt(dictelem[l][d]['bgal'], dictelem[l][d]['lgal'])
                                except:
                                    dictelem[l][d]['sbrt0018'] = dictelem[l][d]['bgal'] * 0.

                    if elemtype[l] == 'lens':
                        #### distance to the source
                        if lensmodltype != 'none':
                            dictelem[l][d]['diss'] = retr_angldist(gdat, dictelem[l][d]['lgal'],  dictelem[l][d]['bgal'], lgalsour[d], bgalsour[d])
                        
                        if lensmodltype == 'elem' or lensmodltype == 'full':
                            dictelem[l][d]['deflprof'] = empty((gdat.numbanglfull, numbelem[l][d]))
                            dictelem[l][d]['mcut'] = empty(numbelem[l][d])
                            dictelem[l][d]['rele'] = empty(numbelem[l][d])
                            dictelem[l][d]['reln'] = empty(numbelem[l][d])
                            dictelem[l][d]['reld'] = empty(numbelem[l][d])
                            dictelem[l][d]['relc'] = empty(numbelem[l][d])
                            for k in arange(numbelem[l][d]):
                                
                                if gdat.variasca:
                                    asca = dictelem[l][d]['asca'][k]
                                else:
                                    asca = gdat.ascaglob

                                if gdat.variacut:
                                    acut = dictelem[l][d]['acut'][k]
                                else:
                                    acut = gdat.acutglob
                                
                                #### deflection profiles
                                dictelem[l][d]['deflprof'][:, k] = retr_deflcutf(gdat.meananglfull, dictelem[l][d]['defs'][k], asca, acut)
             
                                ### truncated mass 
                                dictelem[l][d]['mcut'][k] = retr_mcut(gdat, dictelem[l][d]['defs'][k], asca, acut)
                                
                                #### dot product with the source flux gradient
                                # temp -- weigh the energy and PSF bins
                                dictelem[l][d]['rele'][k] = retr_rele(gdat, cntp['lens'][d, 0, :, 0], dictelem[l][d]['lgal'][k], dictelem[l][d]['bgal'][k], \
                                                                                                                                        1., asca, acut, gdat.indxpixl)
                                dictelem[l][d]['reln'][k] = dictelem[l][d]['rele'][k] / dictelem[l][d]['defs'][k] * gdat.sizepixl
                                dictelem[l][d]['reld'][k] = retr_rele(gdat, gdat.cntpdata[d, 0, :, 0], dictelem[l][d]['lgal'][k], dictelem[l][d]['bgal'][k], 1., \
                                                                                                                                    asca, acut, gdat.indxpixl)
                                dictelem[l][d]['relc'][k] = -retr_rele(gdat, cntp['lens'][d, 0, :, 0], dictelem[l][d]['lgal'][k], dictelem[l][d]['bgal'][k], 1., \
                                                                                                                                    asca, acut, gdat.indxpixl, absv=False)
                   
            ### distribution of element parameters and features
            #### calculate the model filter
            listindxelemfilt = [[[[] for d in indxregipopl[l]] for l in indxpopl] for namefilt in gdat.listnamefilt]
            for k, namefilt in enumerate(gdat.listnamefilt):
                for l in indxpopl:
                    for d in indxregipopl[l]:
                        if namefilt == '':
                            listindxelemfilt[k][l][d] = arange(numbelem[l][d])
                        if namefilt == 'imagbndr':
                            listindxelemfilt[k][l][d] = where((fabs(dictelem[l][d]['lgal']) < gdat.maxmgangdata) & (fabs(dictelem[l][d]['bgal']) < gdat.maxmgangdata))[0]
                        if namefilt == 'deltllik':
                            listindxelemfilt[k][l][d] = where(dictelem[l][d]['deltllik'] > 0.5 * numbcomp[l])[0]
                        if namefilt == 'nrel':
                            listindxelemfilt[k][l][d] = where(dictelem[l][d]['reln'] > 0.3)[0]
    
            listnamefeatsele = getattr(gdat, strgmodl + 'listnamefeatsele')
            
            # histograms of element features
            for l in indxpopl:
                for d in indxregipopl[l]:
                    #### one dimensional
                    for strgfeat in liststrgfeat[l]:
                        if strgfeat == 'spec':
                            temp = zeros((gdat.numbbinsplot, gdat.numbener))
                        else:
                            temp = zeros(gdat.numbbinsplot)
                        dictelem[l][d]['hist' + strgfeat] = temp
                        if strgfeat == 'specplot' or strgfeat == 'deflprof':
                            continue
                        elif strgfeat == 'spec':
                            for i in gdat.indxener:
                                dictelem[l][d]['hist' + strgfeat][:, i] = histogram(dictelem[l][d]['spec'][i, listindxelemfilt[0][l][d]], gdat.binsspec)[0]
                        elif strgfeat == 'cnts':
                            dictelem[l][d]['hist' + strgfeat] = histogram(dictelem[l][d]['cnts'][listindxelemfilt[0][l][d]], gdat.binscnts)[0]
                        elif not (strgfeat == 'curv' and spectype[l] != 'curv' or strgfeat == 'expc' and spectype[l] != 'expc' or strgfeat.startswith('sindarry') and \
                                                                                                                                                        spectype[l] != 'colr'):
                            # temp -- implement the other filters 
                            bins = getattr(gdat, 'bins' + strgfeat)
                            if strgfeat[:-4] in listnamefeatsele[l] and strgmodl == 'true':
                                # temp
                                pass
                                #indx = intersect1d(listindxelemfilt[1][l][d], listindxelemfilt[0][l][d])
                                #if indx.size > 0:
                                #    dictelem[l][d]['hist' + strgfeat] = histogram(dictelem[l][d][strgfeat[:-4]][indx], bins)[0]
                            else:
                                if len(dictelem[l][d][strgfeat]) > 0 and len(listindxelemfilt[0][l][d]) > 0:
                                    # temp
                                    try:
                                        dictelem[l][d]['hist' + strgfeat] = histogram(dictelem[l][d][strgfeat][listindxelemfilt[0][l][d]], bins)[0]
                                    except:
                                        print 'hey'
                                        print 'histograming failed'
                                        print 'strgfeat'
                                        print strgfeat
                                        print 'listindxelemfilt'
                                        print listindxelemfilt
                                        print 'dictelem[l][d][strgfeat]'
                                        print dictelem[l][d][strgfeat]
                                        print 'bins'
                                        print bins
                                        print
                                        raise Exception('')
                
                    #### two dimensional
                    for a, strgfrst in enumerate(liststrgfeatcorr[l]):
                        for b, strgseco in enumerate(liststrgfeatcorr[l]):
                            dictelem[l][d]['hist' + strgfrst + strgseco] = zeros((gdat.numbbinsplot, gdat.numbbinsplot))
                            if a < b:
                                binsfrst = getattr(gdat, 'bins' + strgfrst)
                                binsseco = getattr(gdat, 'bins' + strgseco)
                                if len(dictelem[l][d][strgfrst]) > 0 and len(dictelem[l][d][strgseco]) > 0:
                                    dictelem[l][d]['hist' + strgfrst + strgseco] = histogram2d(dictelem[l][d][strgfrst][listindxelemfilt[0][l][d]], \
                                                                                                 dictelem[l][d][strgseco][listindxelemfilt[0][l][d]], [binsfrst, binsseco])[0]
                                setattr(gdatobjt, strgpfix + 'hist' + strgfrst + strgseco + 'pop%dreg%d' % (l, d), dictelem[l][d]['hist' + strgfrst + strgseco])
                
                    ### priors on element parameters and features
                    dictelem[l][d]['hist' + strgfeat + 'prio'] = empty(gdat.numbbinsplotprio)
                    for strgfeat in liststrgfeatprio[l]:
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
                            if strgfeat == 'expc' and spectype[l] == 'expc':
                                pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                            else:
                                pdfn = pdfn_gaus(xdat, meanvarb, stdv)
                            booltemp = True
                        
                        if booltemp:
                            dictelem[l][d]['hist' + strgfeat + 'prio'] = meanelem[l] * pdfn * deltprio * delt[0] / deltprio[0]
                        
                        setattr(gdatobjt, strgpfix + 'hist' + strgfeat + 'pop%dreg%dprio' % (l, d), dictelem[l][d]['hist' + strgfeat + 'prio'])
                        if strgmodl == 'true':
                            setattr(gdatobjt, 'refrhist' + strgfeat + 'pop%dreg%dprio' % (l, d), dictelem[l][d]['hist' + strgfeat + 'prio'])
        
        if numbtrap > 0:
            for l in indxpopl:
                if elemtype[l] == 'lens':
                    if numbelempopl[l] > 0:
                        ## total truncated mass of the subhalo as a cross check
                        # temp -- generalize
                        for d in indxregipopl[l]:
                            if gdat.variasca:
                                asca = dictelem[l][d]['asca']
                            else:
                                asca = gdat.ascaglob
                            if gdat.variacut:
                                acut = dictelem[l][d]['acut']
                            else:
                                acut = gdat.acutglob
                            factmcutfromdefs = retr_factmcutfromdefs(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour, asca, acut) 
                            masssubh = array([sum(factmcutfromdefs * dictelem[l][d]['defs'])])
            
        ## calculate the host mass, subhalo mass and its fraction as a function of host halo-centric angle and interpolate at the Einstein radius of the host
        if lensmodltype != 'none' or numbtrap > 0 and 'lens' in elemtype:
            listnametemp = ['delt', 'intg']
            listnamevarbmass = []
            listnamevarbmassscal = []
            listnamevarbmassvect = []
            
        if lensmodltype == 'host' or lensmodltype == 'full':
            listnamevarbmassscal += ['masshosttotl']
            for strgtemp in listnametemp:
                listnamevarbmassvect.append('masshost' + strgtemp)
                listnamevarbmassscal.append('masshost' + strgtemp + 'bein')
        if numbtrap > 0 and 'lens' in elemtype:
            listnamevarbmassscal.append('masssubhtotl')
            listnamevarbmassscal.append('fracsubhtotl')
            for strgtemp in listnametemp:
                listnamevarbmassvect.append('masssubh' + strgtemp)
                listnamevarbmassvect.append('fracsubh' + strgtemp)
                listnamevarbmassscal.append('masssubh' + strgtemp + 'bein')
                listnamevarbmassscal.append('fracsubh' + strgtemp + 'bein')

        if lensmodltype != 'none' or numbtrap > 0 and 'lens' in elemtype:
            # find the host, subhalo masses and subhalo mass fraction as a function of halo-centric radius
            for d in gdat.indxregi:
                strgregi = 'reg%d' % d
                angl = sqrt((gdat.meanlgalcart - lgalhost[d])**2 + (gdat.meanbgalcart - bgalhost[d])**2)
            
                for name in listnamevarbmassvect:
                    dicttert[d][name] = zeros(gdat.numbanglhalf)
                    for k in gdat.indxanglhalf:
                        if name[4:8] == 'host':
                            convtemp = conv[d, :]
                        if name[4:8] == 'subh':
                            convtemp = convelem[d, :]
                        
                        if name.endswith('delt'):
                            indxpixl = where((gdat.binsanglhalf[k] < angl) & (angl < gdat.binsanglhalf[k+1]))
                            dicttert[d][name][k] = sum(convtemp[indxpixl]) * gdat.mdencrit * gdat.apix * gdat.adishost**2# / 2. / pi * gdat.deltanglhalf[k]
                        if name.endswith('intg'):
                            indxpixl = where(angl < gdat.meananglhalf[k])
                            dicttert[d][name][k] = sum(convtemp[indxpixl]) * gdat.mdencrit * gdat.apix * gdat.adishost**2# * indxpixl[0].size
                        
                        if name[:4] == 'frac':
                            if dicttert[d]['masshost' + name[8:]][k] != 0.:
                                dicttert[d]['fracsubh' + name[8:]][k] = dicttert[d]['masssubh' + name[8:]][k] / dicttert[d]['masshost' + name[8:]][k]
                    setattr(gdatobjt, strgpfix + name + strgregi, dicttert[d][name])
                    
                    if False:
                        print 'name'
                        print name
                        print 'getattr(gdatobjt, strgpfix + name + strgregi)'
                        print getattr(gdatobjt, strgpfix + name + strgregi)
                        print
            
                # interpolate the host, subhalo masses and subhalo mass fraction at the Einstein radius and save it as a scalar variable
                for name in listnamevarbmassvect:
                    dicttert[d][name + 'bein'] = interp(beinhost[d], gdat.meananglhalf, dicttert[d][name])
                    setattr(gdatobjt, strgpfix + name + 'bein' + strgregi, array([dicttert[d][name + 'bein']]))
                    
                    if False:
                        print 'name'
                        print name
                        print 'dicttert[d][name]'
                        summgene(dicttert[d][name])
                        print 'getattr(gdatobjt, strgpfix + name + bein + strgregi)'
                        print getattr(gdatobjt, strgpfix + name + 'bein' + strgregi)
                        print 
            
        if numbtrap > 0:
            ## copy element features to the global object
            feat = [[[] for d in indxregipopl[l]] for l in indxpopl]
            for l in indxpopl:
                for d in indxregipopl[l]:
                    feat[l][d] = dict()
                    for strgfeat in liststrgfeat[l]:
                        if len(dictelem[l][d][strgfeat]) > 0:
                            if strgmodl == 'true':
                                shap = list(ones(dictelem[l][d][strgfeat].ndim, dtype=int))
                                feat[l][d][strgfeat] = tile(dictelem[l][d][strgfeat], [3] + shap)
                            if strgmodl == 'fitt':
                                feat[l][d][strgfeat] = dictelem[l][d][strgfeat]
                        setattr(gdatobjt, strgpfix + 'hist' + strgfeat + 'pop%dreg%d' % (l, d), dictelem[l][d]['hist' + strgfeat])
                        
            for strgfeat in liststrgfeattotl:
                feattemp = [[[] for d in indxregipopl[l]] for l in indxpopl]
                for l in indxpopl:
                    for d in indxregipopl[l]:
                        if strgfeat in liststrgfeat[l]:
                            if strgfeat in feat[l][d]:
                                feattemp[l][d] = feat[l][d][strgfeat]
                            else:
                                feattemp[l][d] = array([])
                setattr(gdatobjt, strgpfix + strgfeat, feattemp)
            
            if strgmodl == 'true':
                for q in gdat.indxrefr:
                    if elemtype[q] == 'lens' and numbelempopl[q] > 0:
                        for d in gdat.trueindxregipopl[q]:
                            refrmcut = dictelem[q][d]['mcut']
                            refrhistmcut = getattr(gdat, 'truehistmcutpop%dreg%d' % (q, d))
                            indx = where(refrhistmcut > 0)[0]
                            histmcutcorr = zeros(gdat.numbbinsplot)
                            print 'qd'
                            print q, d
                            print 'listindxelemfilt'
                            print listindxelemfilt
                            print 'refrmcut'
                            print refrmcut
                            if len(refrmcut) > 0:
                                refrhistmcutpars = histogram(refrmcut[listindxelemfilt[0]], bins=bins)[0]
                                if indx.size > 0:
                                    histmcutcorr[indx] = refrhistmcutpars[indx] / refrhistmcut[indx]
                            setattr(gdatobjt, strgpfix + 'histmcutcorr', histmcutcorr)
        
        # copy true state to the reference state
        if strgmodl == 'true':
            for name, valu in deepcopy(gdat.__dict__).iteritems():
                if name.startswith('true'):
                        
                    indx = name.find('pop')
                    if indx != -1 and not name.endswith('pop') and name[indx+3].isdigit():
                        if name.startswith('truehist'):
                            print 'replacing'
                        namerefr = name.replace('pop%s' % name[indx+3], 'ref%s' % name[indx+3])
                    else:
                        namerefr = name
                    
                    namerefr = namerefr.replace('true', 'refr')
                    setattr(gdat, namerefr, valu)
        
        ### Exculusive comparison with the true state
        if strgmodl != 'true' and gdat.datatype == 'mock':
             
            if 'lens' in elemtype:
                numbsingcomm = min(deflsing.shape[3], gdat.truedeflsing.shape[3])
                setattr(gdatobjt, strgpfix + 'numbsingcomm', numbsingcomm)
                deflsingresi = deflsing[..., :numbsingcomm] - gdat.truedeflsing[..., :numbsingcomm]
                deflresi = defl - gdat.truedefl
                deflresimgtd = sqrt(sum(deflresi**2, axis=2))
                deflresiperc = 100. * deflresimgtd / gdat.truedeflmgtd
                deflsingresimgtd = sqrt(sum(deflsingresi**2, axis=2))
                deflsingresiperc = 100. * deflsingresimgtd / gdat.truedeflsingmgtd[..., :numbsingcomm]
                setattr(gdatobjt, strgpfix + 'deflsingresi', deflsingresi)
                setattr(gdatobjt, strgpfix + 'deflresi', deflresi)
                setattr(gdatobjt, strgpfix + 'deflresimgtd', deflresimgtd)
                if numbtrap > 0:
                    convelemresi = convelem - gdat.trueconvelem
                    convelemresiperc = 100. * convelemresi / gdat.trueconvelem
                    setattr(gdatobjt, strgpfix + 'convelemresi', convelemresi)
                    setattr(gdatobjt, strgpfix + 'convelemresiperc', convelemresiperc)
                magnresi = magn - gdat.truemagn
                magnresiperc = 100. * magnresi / gdat.truemagn
                setattr(gdatobjt, strgpfix + 'magnresi', magnresi)
                setattr(gdatobjt, strgpfix + 'magnresiperc', magnresiperc)
    
        if numbtrap > 0:
            if gdat.allwrefr:
                # correlate the catalog sample with the reference catalog
                if gdat.refrinfo and not (strgmodl == 'true' and gdat.datatype == 'mock') and gdat.asscrefr:
                    for q in gdat.indxrefr:
                        for l in gdat.fittindxpopl:
                            for d in gdat.fittindxregipopl[l]:
                                if gdat.refrnumbelem[q][d] > 0:
                                    cmpl = array([float(len(indxelemrefrasschits[q][l][d])) / gdat.refrnumbelem[q][d]])
                                    if gdat.diagmode:
                                        if cmpl > 1. or cmpl < 0.:
                                            raise Exception('')
                                else:
                                    cmpl = array([-1.])
                                setattr(gdatobjt, strgpfix + 'cmplref%dpop%dreg%d' % (q, l, d), cmpl)
                                if numbelem[l][d] > 0:
                                    fdis = array([float(indxelemfittasscfals[q][l][d].size) / numbelem[l][d]])
                                    if gdat.diagmode:
                                        if fdis > 1. or fdis < 0.:
                                            raise Exception('')
                                else:
                                    fdis = array([-1.])
                                setattr(gdatobjt, strgpfix + 'fdisref%dpop%dreg%d' % (q, l, d), fdis)
                                
                    # collect the associated fitting element feature for each reference element
                    featrefrassc = [[[[] for d in gdat.fittindxregipopl[l]] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                    for q in gdat.indxrefr:
                        for l in gdat.fittindxpopl:
                            for d in gdat.fittindxregipopl[l]:
                                featrefrassc[q][l][d] = dict()
                                for strgfeat in gdat.refrliststrgfeat[q]:
                                    if not strgfeat in liststrgfeat[l] or strgfeat in gdat.refrliststrgfeatonly[q][l]:
                                        continue
                                    if isinstance(dictelem[l][d][strgfeat], ndarray) and dictelem[l][d][strgfeat].ndim > 1:
                                        continue
                                    featrefrassc[q][l][d][strgfeat] = zeros(gdat.refrnumbelem[q][d]) 
                                    if len(indxelemrefrasschits[q][l][d]) > 0 and len(dictelem[l][d][strgfeat]) > 0:
                                        featrefrassc[q][l][d][strgfeat][indxelemrefrasschits[q][l][d]] = dictelem[l][d][strgfeat][indxelemfittasschits[q][l][d]]
                                    name = strgpfix + strgfeat + 'asscref%dpop%dreg%d' % (q, l, d)
                                    setattr(gdatobjt, name, featrefrassc[q][l][d][strgfeat])

                    # completeness
                    for q in gdat.indxrefr:
                        for l in gdat.fittindxpopl:
                            for d in gdat.fittindxregipopl[l]:
                                if gdat.refrnumbelem[q][d] == 0:
                                    continue
                                for strgfeat in gdat.refrliststrgfeat[q]:
                                    cmplfeat = zeros(gdat.numbbinsplot) - 1.
                                    errrcmplfeat = zeros((2, gdat.numbbinsplot))
                                    refrfeat = getattr(gdat, 'refr' + strgfeat)
                                    if not strgfeat in liststrgfeatodim[l]:
                                        continue
                                    refrhistfeat = getattr(gdat, 'refrhist' + strgfeat + 'ref%dreg%d' % (q, d))
                                    bins = getattr(gdat, 'bins' + strgfeat)
                                    if len(indxelemrefrasschits[q][l][d]) > 0:
                                        histfeatrefrassc = histogram(refrfeat[q][d][0, indxelemrefrasschits[q][l][d]], bins=bins)[0]
                                        indxgood = where(refrhistfeat != 0.)[0]
                                        if indxgood.size > 0:
                                            cmplfeat[indxgood] = histfeatrefrassc[indxgood].astype(float) / refrhistfeat[indxgood]
                                            errrcmplfeat[:, indxgood] = (cmplfeat[indxgood] / sqrt(maximum(ones(indxgood.size), refrhistfeat[indxgood])))[None, :]
                                            if gdat.diagmode:
                                                if where((cmplfeat[indxgood] > 1.) | (cmplfeat[indxgood] < 0.))[0].size > 0:
                                                    raise Exception('')
                                    
                                    setattr(gdatobjt, strgpfix + 'cmpl' + strgfeat + 'ref%dpop%dreg%d' % (q, l, d), cmplfeat)
                                    setattr(gdatobjt, strgpfix + 'errrcmpl' + strgfeat + 'ref%dpop%dreg%d' % (q, l, d), errrcmplfeat)
                       
                    # false discovery rate
                    for q in gdat.indxrefr:
                        for l in gdat.fittindxpopl:
                            for d in gdat.fittindxregipopl[l]:
                                for strgfeat in liststrgfeatodim[l]:
                                    fdisfeat = zeros(gdat.numbbinsplot) - 1.
                                    errrfdisfeat = zeros((2, gdat.numbbinsplot))
                                    if not strgfeat in liststrgfeatodim[l] or strgfeat in gdat.refrliststrgfeatonly[q][l]:
                                        continue
                                    bins = getattr(gdat, 'bins' + strgfeat)
                                    if len(indxelemfittasscfals[q][l][d]) > 0 and len(dictelem[l][d][strgfeat]) > 0:
                                        histfeatfals = histogram(dictelem[l][d][strgfeat][indxelemfittasscfals[q][l][d]], bins=bins)[0]
                                        fitthistfeat = getattr(gdatobjt, strgpfix + 'hist' + strgfeat + 'pop%dreg%d' % (l, d))
                                        indxgood = where(fitthistfeat != 0.)[0]
                                        if indxgood.size > 0:
                                            fdisfeat[indxgood] = histfeatfals[indxgood].astype(float) / fitthistfeat[indxgood]
                                            errrfdisfeat[:, indxgood] = (fdisfeat[indxgood] / sqrt(maximum(ones(indxgood.size), fitthistfeat[indxgood])))[None, :]
                                            if gdat.diagmode:
                                                if where((fdisfeat[indxgood] > 1.) | (fdisfeat[indxgood] < 0.))[0].size > 0:
                                                    raise Exception('')
                                    setattr(gdatobjt, strgpfix + 'fdis' + strgfeat + 'ref%dpop%dreg%d' % (q, l, d), fdisfeat)
                                    setattr(gdatobjt, strgpfix + 'errrfdis' + strgfeat + 'ref%dpop%dreg%d' % (q, l, d), errrfdisfeat)
    
    
            namefeatampl = getattr(gdat, strgmodl + 'namefeatampl')
          
            # temp
            if strgmodl == 'true' and gdat.verbtype > 0:
                for l in indxpopl:
                    for d in gdat.trueindxregipopl[l]:
                        for strgfeat in liststrgfeat[l]:
                            minm = getattr(gdat, 'minm' + strgfeat)
                            maxm = getattr(gdat, 'maxm' + strgfeat)
                            if where(minm > dictelem[l][d][strgfeat])[0].size > 0 or where(maxm < dictelem[l][d][strgfeat])[0].size > 0:
                                print 'Warning: element feature outside the plot limits.'
                                print 'Feature: '
                                print strgfeat
                                print 'Plot minmimum'
                                print minm
                                print 'Plot maxmimum'
                                print maxm
                                if len(dictelem[l][d][strgfeat]) > 0:
                                    print 'Feature minimum'
                                    print amin(dictelem[l][d][strgfeat])
                                    print 'Feature maximum'
                                    print amax(dictelem[l][d][strgfeat])
                                    print
                                if strgfeat == namefeatampl[l]:
                                    raise Exception('')


def proc_cntpdata(gdat):

    # exclude voxels with vanishing exposure
    ## data counts
    if gdat.datatype == 'inpt':
        gdat.cntpdata = retr_cntp(gdat, gdat.sbrtdata, gdat.indxregi, gdat.indxcubeeval)
   
    # correct the likelihoods for the constant data dependent factorial
    gdat.llikoffs = [[] for d in gdat.indxregi]
    for d in gdat.indxregi:
        gdat.llikoffs[d] = sum(sp.special.gammaln(gdat.cntpdata[d] + 1))
    
    if not gdat.killexpo and amax(gdat.cntpdata) == 0:
        print 'gdat.deltener'
        summgene(gdat.deltener)
        print 'gdat.expo'
        summgene(gdat.expo)
        print 'gdat.cntpdata'
        summgene(gdat.cntpdata)
        raise Exception('Data is zero.')

    ## spatial average
    gdat.sbrtdatamean = retr_spatmean(gdat, gdat.cntpdata, gdat.indxregi, boolcntp=True)
    
    # temp
    #gdat.sbrtdatabrod = sum(sum(gdat.cntpdata, axis=2), axis=0)
    
    retr_datatick(gdat)
    
    # 1-point function of the data counts
    for d in gdat.indxregi: 
        for m in gdat.indxevtt:
            if gdat.numbpixl > 1:
                for i in gdat.indxener: 
                    histcntp = histogram(gdat.cntpdata[d][i, :, m], bins=gdat.binscntpdata)[0]
                    setattr(gdat, 'histcntpdatareg%dene%devt%d' % (d, i, m), histcntp)
            else:
                histcntp = histogram(gdat.cntpdata[d][:, 0, m], bins=gdat.binscntpdata)[0]
                setattr(gdat, 'histcntpdatareg%devt%d' % (d, m), histcntp)

    # obtain cartesian versions of the maps
    if gdat.pixltype == 'cart':
        ## data counts
        cntpdatacarttemp = zeros((gdat.numbregi, gdat.numbener, gdat.numbpixlfull, gdat.numbevtt))
        cntpdatacarttemp[:, :, gdat.indxpixlrofi, :] = gdat.cntpdata
        gdat.cntpdatacart = cntpdatacarttemp.reshape((gdat.numbregi, gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
   

def retr_info(pdfnpost, pdfnprio):
    
    info = pdfnpost * log(pdfnpost / pdfnprio)

    return info


def retr_llik_depr(gdat, cntpmodl, indxregieval, indxcubeeval):
   
    llik = [[] for d in indxregieval]
    for d, dd in enumerate(indxregieval):
        llik[dd] = gdat.cntpdata[d][indxcubeeval[0][dd]] * log(cntpmodl[dd]) - cntpmodl[dd]
        if gdat.liketype == 'pois':
            if gdat.verbtype > 1:
                print 'retr_llik_depr'
                print 'gdat.cntpdata[d]'
                summgene(gdat.cntpdata[d])
                print 'cntpmodl[dd]'
                summgene(cntpmodl[dd])
                print 'llik[dd]'
                summgene(llik[dd])
                print
    
        if gdat.liketype == 'gaus':
            llik[dd] = -0.5 * (gdat.cntpdata[d][indxcubeeval[0][dd]] - cntpmodl[dd])**2 / gdat.cntpdata[d][indxcubeeval[0][dd]]
    
    return llik


def retr_mapsgaus(gdat, lgal, bgal, spec, size, ellp, angl):
    
    rttrmatr = array([[cos(angl), -sin(angl)], [sin(angl), cos(angl)]])
    icovmatr = array([[1. / ((1. - ellp) * size)**2, 0.], [0., 1. / size**2]])

    posi = array([lgalgrid - lgal, bgalgrid - bgal])
    mapsgaus = flux * exp(-0.5 * sum(posi * tensordot(self.icovmatr, posi, (1,0)), 0)) / size**2 / (1. - ellp)
        
    return mapsgaus


def retr_sbrtsers(gdat, lgalgrid, bgalgrid, lgal, bgal, spec, size, ellp, angl, seri=4.):
    
    lgalrttr = (1. - ellp) * (cos(angl) * (lgalgrid - lgal) - sin(angl) * (bgalgrid - bgal))
    bgalrttr = sin(angl) * (lgalgrid - lgal) + cos(angl) * (bgalgrid - bgal) 

    angl = sqrt(lgalrttr**2 + bgalrttr**2)

    # interpolate pixel-convolved Sersic surface brightness
    if gdat.serstype == 'intp':
        shapinpt = angl.shape 
        inpt = empty(list(shapinpt) + [3])
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
        sbrthalf = 1. / 7.2 / pi / halfsers**2
    else:
        sbrthalf= 1. / 2. / pi / exp(factsers) * factsers**(2 * indxsers) / indxsers / sp.special.gamma(2. * indxsers) / halfsers**2
                
    ## surface brightness profile
    sbrtsers = sbrthalf * exp(-factsers * ((angl / halfsers)**(1. / indxsers) - 1.))
    
    return sbrtsers


def copytdgu(varb):
    
    if isinstance(varb, ndarray):
        return copy(varb)
    else:
        return deepcopy(varb)


def prep_gdatmodi(gdat, gdatmodi, gdatobjt, strgstat, strgmodl):
    
    if gdat.verbtype > 1:
        print 'prep_gdatmodi()'
    
    strgpfixthis = retr_strgpfix('this', strgmodl)

    gdatmodi.thischro = zeros(gdat.numbchro)
    gdatmodi.thisstdpscalfact = 1.
    gdatmodi.thistmprfactstdv = 1.
    for namevarbstat in gdat.listnamevarbstat:
        temp = copytdgu(getattr(gdatobjt, strgpfixthis + namevarbstat))
        setattr(gdatmodi, strgpfixthis + namevarbstat, temp)
    

def make_anim(strgsrch):

    pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
    listpathruns = fnmatch.filter(os.listdir(pathimag), strgsrch)
    
    print 'Making animations of frame plots...'
    
    for pathruns in listpathruns:
        for namesampdist in ['prio', 'post']:
            for nameextn in ['', 'assc/', 'histodim/', 'histtdim/', 'scattdim/']:
                
                pathframextn = pathimag + pathruns + '/' + namesampdist + '/fram/' + nameextn
                pathanimextn = pathimag + pathruns + '/' + namesampdist + '/anim/' + nameextn
            
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
                            indxfileanim = arange(indxfilelowr, numbfile)
                        else:
                            continue
                            
                        indxfileanim = choice(indxfileanim, replace=False, size=indxfileanim.size)
                        
                        cmnd = 'convert -delay 20 -density 300 -quality 100 '
                        for n in range(indxfileanim.size):
                            cmnd += '%s%s ' % (pathframextn, listfile[indxfileanim[n]])
    
                        namegiff = '%s%s.gif' % (pathanimextn, name + liststrgextn[k])
                        cmnd += ' ' + namegiff
                        if not os.path.exists(namegiff):
                            print 'Run: %s, pdf: %s' % (pathruns, namesampdist)
                            print 'Making %s animation...' % name
                            print
                            os.system(cmnd)
                        else:
                            pass
    

