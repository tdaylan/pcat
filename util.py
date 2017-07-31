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


def retr_probpois(data, modl):
    
    lprb = data * log(modl) - modl - sp.special.gammaln(data + 1)
    
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


## proposals
### decode the transdimensional element list
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
                if spectype[l] == 'expc':
                    indxsampcomp['expc'][l] = indxsamptemp + 4
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


### update the current state with a Gaussian draw with potential tail
def retr_gaus(gdat, gdatmodi, indxparamodi, stdvpara):
    
    numbparamodi = indxparamodi.size
    stdvtemp = normal(size=numbparamodi) * stdvpara
    if rand() < gdat.fracproprand:
        stdvtemp[choice(indxparamodi)] += randn()
    
    gdatmodi.nextsamp[indxparamodi] = gdatmodi.thissamp[indxparamodi] + stdvtemp
    gdatmodi.nextsamp[indxparamodi] = gdatmodi.nextsamp[indxparamodi] % 1.

       
### draw proposal type
def retr_thisindxprop(gdat, gdatmodi, thisindxpopl=None, brth=False, deth=False):
    
    gdatmodi.thisaccppsfn = True
    gdatmodi.thisaccpprio = True
    
    gdatmodi.evalllikpert = False 

    gdatmodi.thisdeltlpri = 0.

    # initialize the Boolean flag indicating the type of transdimensional proposal
    gdatmodi.propcomp = False
    gdatmodi.propfixp = False
    gdatmodi.proppsfp = False
    gdatmodi.propwith = False
    gdatmodi.propbrth = False
    gdatmodi.propdeth = False
    gdatmodi.propsplt = False
    gdatmodi.propmerg = False
    
    # temp
    gdatmodi.propmeanelem = False
    gdatmodi.propdist = False
    gdatmodi.prophypr = False

    gdatmodi.thisljcbfact = 0.
    gdatmodi.thislpautotl = 0. 
    gdatmodi.thislcomfact = 0.
   
    gdatmodi.prophost = False
    gdatmodi.proppsfnconv = False
    
    # index of the population in which a transdimensional proposal will be made
    if thisindxpopl == None:
        gdatmodi.indxpoplmodi = choice(gdat.fittindxpopl)
    else:
        gdatmodi.indxpoplmodi = thisindxpopl

    if (rand() < gdat.probtran or brth or deth) and (gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] != gdat.fittminmnumbelem[gdatmodi.indxpoplmodi] or \
                                                     gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] != gdat.fittmaxmnumbelem[gdatmodi.indxpoplmodi]):
        
        if rand() < gdat.probbrde or brth or deth:
            
            ## births and deaths
            if gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] == gdat.fittmaxmnumbelem[gdatmodi.indxpoplmodi] or deth:
                gdatmodi.propdeth = True
            elif gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] == gdat.fittminmnumbelem[gdatmodi.indxpoplmodi] or brth:
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
            if gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] == gdat.fittmaxmnumbelem[gdatmodi.indxpoplmodi]:
                gdatmodi.propmerg = True
            elif gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] == gdat.fittminmnumbelem[gdatmodi.indxpoplmodi]:
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
            
            if gdatmodi.indxsampmodi in gdat.fittindxfixp:
                gdatmodi.propfixp = True
            else:
                gdatmodi.propcomp = True
                
            if gdat.fittpsfnevaltype != 'none':
                if gdatmodi.indxsampmodi in gdat.fittindxfixppsfp:
                    gdatmodi.proppsfp = True

            # temp
            if gdat.fittnumbtrap > 0:
                if gdatmodi.indxsampmodi in gdat.fittindxfixpmeanelem:
                    gdatmodi.propmeanelem = True
                if gdatmodi.indxsampmodi in gdat.fittindxfixpdist:
                    gdatmodi.propdist = True

            gdatmodi.thisindxproptype = gdat.indxstdppara[gdatmodi.indxsampmodi]
        else:
            gdatmodi.thisindxproptype = gdat.indxproptypewith
        gdatmodi.propwith = True
        gdatmodi.proppsfnconv = (gdat.fittpsfnevaltype == 'conv' or gdat.fittpsfnevaltype == 'full') and gdatmodi.proppsfp
        if gdat.elemtype == 'lens':
            gdatmodi.prophost = gdatmodi.indxsampmodi in gdat.fittindxfixphost
    
    gdatmodi.proptran = gdatmodi.propbrth or gdatmodi.propdeth or gdatmodi.propsplt or gdatmodi.propmerg
    
    gdatmodi.prophypr = gdatmodi.propmeanelem or gdatmodi.propdist
    
    gdatmodi.evalllikpert = (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and (not gdatmodi.propfixp or gdatmodi.proppsfp and gdat.fittnumbtrap > 0)
    
    if gdat.verbtype > 1:
        print 'retr_thisindxprop():'
        if gdat.fittnumbtrap > 0:
            print 'gdatmodi.thisindxelemfull'
            print gdatmodi.thisindxelemfull
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
        print 'gdatmodi.indxpoplmodi'
        print gdatmodi.indxpoplmodi
        print 'gdatmodi.evalllikpert'
        print gdatmodi.evalllikpert
        print 'gdatmodi.propfixp'
        print gdatmodi.propfixp
        print 'gdatmodi.proppsfp'
        print gdatmodi.proppsfp
        print 'gdatmodi.prophost'
        print gdatmodi.prophost
        print 'gdatmodi.propcomp'
        print gdatmodi.propcomp
        if gdatmodi.propwith and gdat.propwithsing:
            print 'indxsampmodi'
            print gdatmodi.indxsampmodi
        print
    
    if gdat.diagmode:
        if gdat.probbrde == 0. and (gdatmodi.propbrth or gdatmodi.propdeth):
            raise Exception('')


### update sampler state
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
        
    if gdat.penalpridiff and not gdatmodi.prophypr:
        gdatmodi.thislpridiff = gdatmodi.nextlpridiff
        
    if gdat.calcllik:
        # update the log-likelihood
        gdatmodi.thislliktotl = copy(gdatmodi.nextlliktotl)
        
        # update the log-likelihood
        gdatmodi.thisllik[:, gdatmodi.indxpixleval, :] = copy(gdatmodi.nextllik)

        if gdatmodi.proppsfnconv:
            gdatmodi.thispsfnconv = [[[] for i in gdat.indxener] for m in gdat.indxevtt]
            for m in gdat.indxevtt:
                for i in gdat.indxener:
                    gdatmodi.thispsfnconv[m][i] = AiryDisk2DKernel(gdatmodi.nextpsfp[i] / gdat.sizepixl)
    
        if gdat.fitthostemistype != 'none':
            if gdatmodi.prophost:
                gdatmodi.thissbrthost = copy(gdatmodi.nextsbrthost)
        
        if gdat.fittlensmodltype != 'none':
            if gdatmodi.prophost:
                gdatmodi.thisdeflhost = copy(gdatmodi.nextdeflhost)

    # temp
    # keep the state vectors clean
    if gdatmodi.propdeth:
        gdatmodi.thissampvarb[gdatmodi.indxsamptran[0]] = 0.
        gdatmodi.thissamp[gdatmodi.indxsamptran[0]] = 0.
    
    if gdat.fittnumbtrap > 0:
        gdatmodi.thisindxelemfull = deepcopy(gdatmodi.nextindxelemfull)
        gdatmodi.thisindxsampcomp = deepcopy(gdatmodi.nextindxsampcomp)
    
        if gdat.calcllik and gdat.propwithsing:
            
            if gdat.elemtype == 'lght':
                if gdat.verbtype > 1:
                    print 'gdatmodi.nextsbrtpnts'
                    summgene(gdatmodi.nextsbrtpnts)
                    print 'gdatmodi.thissbrtpnts'
                    summgene(gdatmodi.thissbrtpnts)

                gdatmodi.thissbrtpnts[:, gdatmodi.indxpixleval, :] = copy(gdatmodi.nextsbrtpnts[:, gdatmodi.indxpixleval, :])
            if gdat.elemtype == 'lens':
                gdatmodi.thisdeflelem = copy(gdatmodi.nextdeflelem)


def initcompfromrefr(gdat, gdatmodi, namerefr):
    
    for l in gdat.fittindxpopl:
        for k, strgcomp in enumerate(gdat.fittliststrgcomp[l]):
            minm = getattr(gdat, 'fittminm' + strgcomp)
            maxm = getattr(gdat, 'fittmaxm' + strgcomp)
            comp = getattr(gdat, namerefr + strgcomp)[l]
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
                compunit = rand(gdat.truenumbelem[l])
            gdatmodi.thissamp[gdatmodi.thisindxsampcomp[strgcomp][l]] = compunit


## model evaluation

### find the spectra of sources
def retr_spec(gdat, flux, sind, curv=None, expc=None, spectype='powr', plot=False):
    
    if gdat.numbener == 1:
        spec = flux[None, :]
    else:
        if plot:
            meanener = gdat.meanenerplot
        else:
            meanener = gdat.meanener

        if spectype == 'powr':
            spec = flux[None, :] * (meanener / gdat.enerpivt)[:, None]**(-sind[None, :])
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
            
            if False:
                print 'dist'
                summgene(dist * gdat.anglfact)
                print 'indxoaxitemp'
                print indxoaxitemp
                print 
            
            psfntemp = psfnintp[indxoaxitemp](dist)
        else:
            psfntemp = psfnintp(dist)
    if gdat.kernevaltype == 'bspx':
        pass

    # scale by the PS spectrum
    sbrtpnts = spec[:, None, None] * psfntemp
                
    return sbrtpnts


### find the set of pixels in proximity to a position on the map
def retr_indxpixlevalconc(gdat, dicttemp, evalcirc):

    lgal = dicttemp['lgaleval']
    bgal = dicttemp['bgaleval']
    if gdat.elemtype == 'lght':
        varbeval = abs(dicttemp['speceval'][gdat.indxenerpivt[0], :])
    if gdat.elemtype == 'lens':
        varbeval = dicttemp['defseval']
    if gdat.elemtype == 'clus':
        varbeval = dicttemp['nobjeval']
    numbelemeval = lgal.size
    
    if False:
        if not isscalar(lgal):
            lgal = lgal[0]
            bgal = bgal[0]
            ampl = ampl[0]
    
    if evalcirc == 'locl':  
        listindxpixleval = []
        for k in range(numbelemeval):
            indxpixlpnts = retr_indxpixl(gdat, bgal[k], lgal[k])
            indxfluxproxtemp = digitize(varbeval[k], gdat.binsprox) - 1
            indxpixleval = gdat.indxpixlprox[indxfluxproxtemp][indxpixlpnts]
            if isinstance(indxpixleval, int):
                indxpixleval = gdat.indxpixl
            listindxpixleval.append(indxpixleval)
            if gdat.verbtype > 1:
                print 'k'
                print k
                print 'dicttemp[lgaleval][k]'
                print dicttemp['lgaleval'][k]
                print 'dicttemp[bgaleval][k]'
                print dicttemp['bgaleval'][k]
                print 'varbeval'
                print varbeval
                print 'indxpixleval'
                summgene(indxpixleval)
    
    if numbelemeval > 0:
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
def retr_cntp(gdat, sbrt, indxpixlmean=None):

    if indxpixlmean == None:
        indxpixlmean = arange(sbrt.shape[1])

    cntp = sbrt * gdat.expo[:, indxpixlmean, :] * gdat.apix
    if gdat.enerdiff:
        cntp *= gdat.deltener[:, None, None]
        
    return cntp


def retr_elpsfrac(elpsaxis):
    
    distnorm = sum(((listsamp - gdat.elpscntr[None, :]) / elpsaxis[None, :])**2, axis=1)
    indxsampregu = where(distnorm < 1.)[0]
    thissampfrac = indxsampregu.size / gdat.numbsamp
    vari = (thissampfrac / 0.05 - 1.)**2
    
    return vari


## plotting
### construct path for plots
def retr_plotpath(gdat, gdatmodi, strg, strgplot, nameinte=''):
    
    if strg == 'true' or strg == '':
        path = gdat.pathinit + nameinte + strgplot + '.pdf'
    elif strg == 'post' or strg == 'mlik':
        path = gdat.pathplot + gdat.namesampdist + '/finl/' + nameinte + strg + strgplot + '.pdf'
    elif strg == 'this':
        path = gdat.pathplot + gdat.namesampdist + '/fram/' + nameinte + strg + strgplot + '_swep%09d.pdf' % gdatmodi.cntrswep
    
    return path


### determine the marker size
def retr_mrkrsize(gdat, compampl):
    
    minm = getattr(gdat, 'minm' + gdat.namecompampl) 
    maxm = getattr(gdat, 'maxm' + gdat.namecompampl) 
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
    if gdat.numbener == 5:
        gdat.psfpexpr = array([0.424 / gdat.anglfact, 2.75, 0.424 / gdat.anglfact, 2.59, 0.440 / gdat.anglfact, 2.47, 0.457 / gdat.anglfact, 2.45, 0.529 / gdat.anglfact, 3.72])
    if gdat.numbener == 2:
        gdat.psfpexpr = array([0.427 / gdat.anglfact, 2.57, 0.449 / gdat.anglfact, 2.49])

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
    
    parastrg = ['score', 'gcore', 'stail', 'gtail', 'ntail']
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
            fermform[:, m, k] = interp1d_pick(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)
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
    

def retr_chandata(gdat):
    
    gdat.numbrefr += 1
    gdat.indxrefr = arange(gdat.numbrefr)
    
    gdat.dictrefr = []
    for q in gdat.indxrefr:
        gdat.dictrefr.append(dict())
    
    gdat.legdrefr = ['Xue+2011']
    
    gdat.listnamerefr += ['xu11']
    
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


    with open(gdat.pathinpt + 'chancatl.txt', 'r') as thisfile:
        
        pathfile = gdat.pathinpt + 'Xue2011.fits'
        
        hdun = pf.open(pathfile)
        hdun[0].data

        lgalchan = hdun[1].data['_Glon']
        bgalchan = hdun[1].data['_Glat']
        fluxchansoft = hdun[1].data['SFlux']
        fluxchanhard = hdun[1].data['HFlux']
        objttypechan = hdun[1].data['Otype']

        #G_long = [] #deg
        #G_lat = [] #deg
        #id_number = []
        #off_angle = [] # arcmin
        #flux_cnts = [] # for xray band
        #soft_cnts = [] # 0.5-2
        #hard_cnts = [] # 2-8
        #c_offset = [] #angular offset between xray and optical/NIR components in arcse
        #C_mag = [] # full optical mag?
        #W_mag = [] # 600-700 nm
        #GD_mag = [] # 750-1000 nm from GOODS-S z-band
        #G_mag = [] # 750-1000 nm from GEMS z-band
        #M_mag = [] # 2-3 micron
        #S_mag = [] # 3-4 micron
        #flux_erg_full = [] # erg/(s*cm^2)
        #flux_erg_soft = []
        #flux_erg_hard = []
        #Otype = [] # AGN/Galaxy/Star
        #for line in thisfile:
        #    line = line.split()
        #    G_long.append(line[0])
        #    G_lat.append(line[1])
        #    id_number.append(line[2])
        #    off_angle.append(line[3])
        #    flux_cnts.append(line[4])
        #    soft_cnts.append(line[5])
        #    hard_cnts.append(line[6])
        #    c_offset.append(line[7])
        #    C_mag.append(line[8])
        #    W_mag.append(line[9])
        #    GD_mag.append(line[10])
        #    G_mag.append(line[11])
        #    M_mag.append(line[12])
        #    S_mag.append(line[13])
        #    flux_erg_full.append(line[14])
        #    flux_erg_soft.append(line[15])
        #    flux_erg_hard.append(line[16])
        #    Otype.append(line[17])
        #lgalchan = (asarray(G_long)).astype(float)
        #bgalchan = (asarray(G_lat)).astype(float)
        #cntschan = (asarray(flux_cnts)).astype(float)
        #cntschansoft = (asarray(soft_cnts)).astype(float)
        #cntschanhard = (asarray(hard_cnts)).astype(float)
        ##offschan = (asarray(c_offset)).astype(float)
        ##cmagchan = (asarray(C_mag)).astype(float)
        ##wmagchan = (asarray(W_mag)).astype(float)
        ##dmagchan = (asarray(GD_mag)).astype(float)
        ##gmagchan = (asarray(G_mag)).astype(float)
        ##mmagchan = (asarray(M_mag)).astype(float)
        ##smagchan = (asarray(S_mag)).astype(float)
        ##fluxchanfull = (asarray(flux_erg_full)).astype(float)
        #fluxchansoft = (asarray(flux_erg_soft)).astype(float)
        #fluxchanhard = (asarray(flux_erg_hard)).astype(float)
        #objttypechan = (asarray(Otype))
    
    indx = where(objttypechan == 'AGN')[0]
    objttypechan[indx] = 0
    indx = where(objttypechan == 'Galaxy')[0]
    objttypechan[indx] = 1
    indx = where(objttypechan == 'Star')[0]
    objttypechan[indx] = 2
    objttypechan = objttypechan.astype(int)

    if gdat.numbener == 2:
        path = gdat.pathinpt + 'CDFS-4Ms-0p5to2-asca-im-bin1.fits'
    else:
        path = gdat.pathinpt + '0.50-0.91_thresh.img'
    
    # rotate reference elements to the spatial coordinate system of PCAT
    listhdun = ap.io.fits.open(path)
    wcso = ap.wcs.WCS(listhdun[0].header)
    skycobjt = ap.coordinates.SkyCoord("galactic", l=lgalchan, b=bgalchan, unit='deg')
    rascchan = skycobjt.fk5.ra.degree
    declchan = skycobjt.fk5.dec.degree

    if gdat.numbsidecart == 300:
        if gdat.anlytype.startswith('extr'):
            indxpixllgal = 1490
            indxpixlbgal = 1430
        if gdat.anlytype.startswith('home'):
            indxpixllgal = 150
            indxpixlbgal = 150
    else:
        raise Exception('Reference elements cannot be aligned with the spatial axes!')
    
    lgalchan, bgalchan = wcso.wcs_world2pix(rascchan, declchan, 0)
    lgalchan -= indxpixllgal + gdat.numbsidecart / 2
    bgalchan -= indxpixlbgal + gdat.numbsidecart / 2
    lgalchan *= gdat.sizepixl
    bgalchan *= gdat.sizepixl
    
    # this twist is intentional
    gdat.refrbgal = [lgalchan]
    gdat.refrlgal = [bgalchan]
    gdat.refrotyp = [objttypechan]
    
    # temp
    gdat.refrspec = [zeros((3, gdat.numbener, gdat.refrlgal[0].size))]

    for q in gdat.indxrefr:
        gdat.refrlgal[q] = tile(gdat.refrlgal[q], (3, 1)) 
        gdat.refrbgal[q] = tile(gdat.refrbgal[q], (3, 1)) 
        gdat.refrotyp[q] = tile(gdat.refrotyp[q], (3, 1)) 
    
    if gdat.numbener == 2:
        gdat.refrspec[0][0, 0, :] = fluxchansoft * 0.624e9
        gdat.refrspec[0][0, 1, :] = fluxchanhard * 0.624e9 / 16.
    else:
        gdat.refrspec[0][0, :, :] = fluxchansoft[None, :] * 0.624e9

    # temp
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :]
    
    # temp
    gdat.refrspec[0][where(gdat.refrspec[0] < 0.)] = 0.
   
    gdat.refrsind = []
    if gdat.numbener > 1:
        gdat.refrsind.append(-log(gdat.refrspec[0][0, 1, :] / gdat.refrspec[0][0, 0, :]) / log(sqrt(7. / 2.) / sqrt(0.5 * 2.)))
        gdat.refrsind[0][where(logical_not(isfinite(gdat.refrsind[0])))[0]] = 2.
    
    # temp
    if gdat.numbener > 1:
        gdat.refrsind[0] = tile(gdat.refrsind[0], (3, 1)) 


def retr_fermdata(gdat):
    
    gdat.listnamerefr += ['ac15']
    gdat.numbrefr += 1
    
    gdat.legdrefr = ['Acero+2015']

    setattr(gdat, 'lablcurvac15', '%s_{3FGL}' % gdat.lablcurv)
    setattr(gdat, 'lablexpcac15', 'E_{c,3FGL}')
    gdat.indxrefr = arange(gdat.numbrefr)

    gdat.minmcurv = -1.
    gdat.maxmcurv = 1.
    gdat.minmexpc = 0.1
    gdat.maxmexpc = 10.
    for name in gdat.listnamerefr:
        setattr(gdat, 'minmcurv' + name, -1.)
        setattr(gdat, 'maxmcurv' + name, 1.)
        setattr(gdat, 'minmexpc' + name, 0.1)
        setattr(gdat, 'maxmexpc' + name, 10.)
    
    path = gdat.pathdata + 'expr/pnts/gll_psc_v16.fit'
    fgl3 = pf.getdata(path)
   
    gdat.refrlgal = [deg2rad(fgl3['glon'])]
    gdat.refrlgal[0] = ((gdat.refrlgal[0] - pi) % (2. * pi)) - pi
    gdat.refrbgal = [deg2rad(fgl3['glat'])]
    
    gdat.truenumbelemfull = gdat.refrlgal[0].size

    for q in gdat.indxrefr:
        gdat.refrlgal[q] = tile(gdat.refrlgal[q], (3, 1)) 
        gdat.refrbgal[q] = tile(gdat.refrbgal[q], (3, 1)) 
    
    gdat.refrspec = [empty((3, gdat.numbener, gdat.truenumbelemfull))]
    gdat.refrspec[0][0, :, :] = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], \
                                                                                            fgl3['Flux10000_100000']))[gdat.indxenerincl, :] / gdat.deltener[:, None]
    

    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], \
                                                        fgl3['Unc_Flux10000_100000']))[gdat.indxenerincl, :, :] / gdat.deltener[:, None, None] 
    gdat.refrspec[0][1, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 0]
    gdat.refrspec[0][2, :, :] = gdat.refrspec[0][0, :, :] + fgl3specstdvtemp[:, :, 1]
    gdat.refrspec[0][where(isfinite(gdat.refrspec[0]) == False)] = 0.
    
    gdat.refrflux = [gdat.refrspec[0][:, gdat.indxenerpivt[0], :]]
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3strg = fgl3['Source_Name']
    fgl3strgclss = fgl3['CLASS1']
    fgl3strgassc = fgl3['ASSOC1']
    
    fgl3timevari = fgl3['Variability_Index']

    fgl3spectype = fgl3['SpectrumType']
    gdat.refrsind = [fgl3['Spectral_Index']]
    gdat.refrcurv = [fgl3['beta']]
    gdat.refrcurv[0][where(logical_not(isfinite(gdat.refrcurv[0])))] = -10.
    gdat.refrexpc = [fgl3['Cutoff'] * 1e-3]
    
    gdat.refrsind[0] = tile(gdat.refrsind[0], (3, 1)) 
    gdat.refrcurv[0] = tile(gdat.refrcurv[0], (3, 1)) 
    gdat.refrexpc[0] = tile(gdat.refrexpc[0], (3, 1)) 
    
    gdat.refrexpc[0][where(logical_not(isfinite(gdat.refrexpc[0])))] = 0.

    indxtimevari = where(fgl3timevari < 100.)[0]
    
    #indxtimevari = where(gdat.refrspec[0][0, gdat.indxenerpivt[0], :] > gdat.minmflux)[0]
    #gdat.refrbgal = gdat.refrbgal[:, indxtimevari]
    #gdat.refrsind = gdat.refrsind[:, indxtimevari]
    #gdat.refrcurv = gdat.refrcurv[:, indxtimevari]
    #gdat.refrexpc = gdat.refrexpc[:, indxtimevari]
    #gdat.refrspec[0] = gdat.refrspec[0][:, :, indxtimevari]


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
            gdatmodi.indxparaprop = concatenate((gdat.indxfixpprop, concatenate(gdatmodi.thisindxsampcomp['comp'])))
            if k in concatenate(gdatmodi.thisindxsampcomp['lgal']):
                print
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


def retr_prop(gdat, gdatmodi, thisindxelem=None):
 
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
            
            if gdatmodi.propfixp:
                gdatmodi.numbelemeval = 0
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
                gdatmodi.thisstdvpara = stdvstdp[gdat.indxstdppara[gdatmodi.indxsampmodi]]
            else:
                retr_propcompscal(gdat, gdatmodi, stdvstdp, gdatmodi.indxpoplmodi, gdatmodi.thisstrgcomp, indxelemfull=gdatmodi.indxelemfullmodi[0])
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdatmodi.thisstdvpara)
            
            if gdatmodi.propfixp:
                gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_fixp(gdat, 'fitt', gdatmodi.nextsamp[gdatmodi.indxsampmodi], gdatmodi.indxsampmodi)
           
            # temp
            if gdat.fittnumbtrap > 0 and gdatmodi.propdist:
                # temp
                # this should check whether rescaling is necessary
                gdatmodi.thisstrgcomp = gdat.fittnamepara[gdatmodi.indxsampmodi][-16:-12]

                ### rescale the element components due to the hyperparameter step
                rscl_elem(gdat, gdatmodi, gdatmodi.indxsampmodi)
            
            if gdatmodi.propcomp:
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
                        retr_propcompscal(gdat, gdatmodi, stdvstdp, l, strgcomp)
                        retr_gaus(gdat, gdatmodi, gdatmodi.thisindxsampcomp[strgcomp][l], gdatmodi.thisstdvpara)

    else:
        # number of unit sample vector elements to be modified
        numbcompmodi = gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
    
    if gdatmodi.proptran:
        gdatmodi.auxipara = rand(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.indxsamptran = []
    
    if gdatmodi.propbrth or gdatmodi.propsplt:
       
        # find an empty slot in the element list
        for k in range(gdat.fittmaxmnumbelem[gdatmodi.indxpoplmodi]):
            if not k in gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi]:
                break
       
        # sample indices to add the new element
        gdatmodi.indxsamptran.append(retr_indxsamppnts(gdat, gdatmodi.indxpoplmodi, array([k])))
        gdatmodi.indxelemmodi = [k]
        gdatmodi.nextindxelemfull[gdatmodi.indxpoplmodi].append(k)
        gdatmodi.indxelemfullmodi = [gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]].astype(int)]
    if gdatmodi.propbrth:
        
        gdatmodi.numbelemeval = 1
        
        # sample auxiliary variables
        gdatmodi.nextsamp[gdatmodi.indxsamptran[0]] = gdatmodi.auxipara
    
    # death
    if gdatmodi.propdeth:
        
        gdatmodi.numbelemeval = 1
        
        # occupied element index to be killed
        if thisindxelem == None:
            dethindxindxelem = choice(arange(gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]], dtype=int))
        else:
            dethindxindxelem = thisindxelem

        # element index to be killed
        gdatmodi.indxelemmodi = []
        gdatmodi.indxelemfullmodi = []
        if gdat.verbtype > 1:
            print 'dethindxindxelem'
            print dethindxindxelem
        gdatmodi.indxelemmodi.append(gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][dethindxindxelem])
        gdatmodi.indxelemfullmodi.append(dethindxindxelem)
        # parameter indices to be killed
        gdatmodi.indxsamptran.append(gdat.fittindxsampcompinit + gdat.fittnumbtrapcuml[gdatmodi.indxpoplmodi] + \
                                                    gdatmodi.indxelemmodi[0] * gdat.fittnumbcomp[gdatmodi.indxpoplmodi] + gdat.fittindxcomp[gdatmodi.indxpoplmodi])
    
    if gdatmodi.propsplt or gdatmodi.propmerg:
        gdatmodi.comppare = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.compfrst = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
        gdatmodi.compseco = empty(gdat.fittnumbcomp[gdatmodi.indxpoplmodi])
    
    # split
    if gdatmodi.propsplt:
        
        gdatmodi.numbelemeval = 3
        
        gdatmodi.indxelemfullsplt = choice(arange(gdatmodi.thissamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]], dtype=int))
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
            gdatmodi.mergindxelemfrst = gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmergfrst]

            ## second element index to be merged
            gdatmodi.mergindxelemseco = gdatmodi.thisindxelemfull[gdatmodi.indxpoplmodi][gdatmodi.indxelemfullmergseco]
            
            # determine indices of the modified elements in the sample vector
            ## first element
            # temp -- this would not work for multiple populations !
            gdatmodi.indxsampfrst = gdat.fittindxsampcompinit + gdat.fittnumbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxelemfrst
            gdatmodi.indxsampfrstfinl = gdatmodi.indxsampfrst + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
            gdatmodi.indxsamptran.append(arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl))

            ## second element
            gdatmodi.indxsampseco = gdat.fittindxsampcompinit + gdat.fittnumbcomp[gdatmodi.indxpoplmodi] * gdatmodi.mergindxelemseco
            gdatmodi.indxsampsecofinl = gdatmodi.indxsampseco + gdat.fittnumbcomp[gdatmodi.indxpoplmodi]
            gdatmodi.indxsamptran.append(arange(gdatmodi.indxsampseco, gdatmodi.indxsampsecofinl))
            
            # indices of the sample vector elements to be modified
            gdatmodi.indxsampmodi = arange(gdatmodi.indxsampfrst, gdatmodi.indxsampfrstfinl)

            # indices of the element to be merged
            gdatmodi.indxelemmodi = [gdatmodi.mergindxelemfrst, gdatmodi.mergindxelemseco]

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
                print 'mergindxfrst: ', gdatmodi.mergindxelemfrst
                print 'gdatmodi.indxelemfullmergfrst: ', gdatmodi.indxelemfullmergfrst
                print 'mergindxseco: ', gdatmodi.mergindxelemseco
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
        gdatmodi.nextsamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] + 1
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] + 1
    if gdatmodi.propdeth or gdatmodi.propmerg:
        gdatmodi.nextsamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] = gdatmodi.thissamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] - 1
        gdatmodi.nextsampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]] - 1
    
    if (gdatmodi.propdeth or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        # remove the element from the occupied element list
        for a, indxelem in enumerate(gdatmodi.indxelemmodi):
            if a == 0 and gdatmodi.propdeth or a == 1 and gdatmodi.propmerg:
                gdatmodi.nextindxelemfull[gdatmodi.indxpoplmodi].remove(indxelem)
    
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
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[0]))
            if gdatmodi.propdeth:
                gdatmodi.indxsampmodi = gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi, None]
            if gdatmodi.propsplt:
                gdatmodi.indxsampmodi = concatenate((gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi, None], gdatmodi.indxsamptran[0], gdatmodi.indxsamptran[1]))
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
        print 'gdatmodi.thisaccpprio'
        print gdatmodi.thisaccpprio
        print 'gdatmodi.indxsampchec'
        print gdatmodi.indxsampchec
        if gdat.fittnumbtrap > 0:
            if gdatmodi.proptran:
                print 'gdatmodi.indxsamptran'
                for indxsamptran in gdatmodi.indxsamptran:
                    print indxsamptran
            print 'gdatmodi.thisindxelemfull'
            print gdatmodi.thisindxelemfull
            print 'gdatmodi.nextindxelemfull'
            print gdatmodi.nextindxelemfull
            if (gdatmodi.propcomp or gdatmodi.proptran) and not (gdatmodi.propmerg and not gdatmodi.thisaccpprio):
                print 'gdatmodi.indxpoplmodi'
                print gdatmodi.indxpoplmodi
                print 'gdatmodi.indxelemmodi'
                print gdatmodi.indxelemmodi
                print 'gdatmodi.indxelemfullmodi'
                print gdatmodi.indxelemfullmodi
    
    if gdat.fittnumbtrap > 0:
        if gdatmodi.propwith:
            gdatmodi.nextindxsampcomp = gdatmodi.thisindxsampcomp
        else:
            gdatmodi.nextindxsampcomp = retr_indxsampcomp(gdat, gdatmodi.nextindxelemfull, 'fitt')
    
    # temp -- this is inefficient for propwithsing proposals
    if gdatmodi.thisaccpprio and (gdatmodi.propbrth or gdatmodi.propsplt or gdatmodi.propwith and not gdat.propwithsing):
        gdatmodi.nextsampvarb = retr_sampvarb(gdat, 'fitt', gdatmodi.nextsamp, gdatmodi.nextindxsampcomp)
   
    if gdat.propwithsing and (gdatmodi.propcomp or gdatmodi.proptran):
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
        
    # calculate the factor, i.e., Jacobian and combinatorial, to multiply the acceptance rate
    if (gdatmodi.propsplt or gdatmodi.propmerg) and gdatmodi.thisaccpprio:
        
        ## Jacobian
        if gdatmodi.propsplt:
            gdatmodi.thisljcbfact = log(gdatmodi.comppare[2])
        else:
            gdatmodi.thisljcbfact = -log(gdatmodi.comppare[2])
         
        ## combinatorial factor
        # temp
        #thisnumbelem = gdatmodi.thissamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]]
        #if gdatmodi.propsplt:
        #    gdatmodi.thislcomfact = log(gdatmodi.thissamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]]**2 / gdatmodi.thisnumbpair)
        #else:
        #    gdatmodi.thislcomfact = log(gdatmodi.thisnumbpair / gdatmodi.thissamp[gdat.fittindxfixpnumbelem[gdatmodi.indxpoplmodi]]**2)
        gdatmodi.thislcomfact = 0.


def retr_propcompscal(gdat, gdatmodi, stdvstdp, l, strgcomp, indxelemfull=slice(None)):
    
    stdvstdpcomp = stdvstdp[getattr(gdat, 'indxstdp' + strgcomp)]
    if False:
        thiscompampl = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
        nextcompampl = gdatmodi.nextsampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
        minmcompampl = getattr(gdat, 'minm' + gdat.namecompampl)
        thiscompunit = gdatmodi.thissamp[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
        nextcompunit = gdatmodi.nextsamp[gdatmodi.thisindxsampcomp[gdat.namecompampl][l][indxelemfull]]
        if strgcomp == gdat.namecompampl:
            # temp -- this only works if compampl is powr distributed
            gdatmodi.thisstdvpara = stdvstdpcomp / (thiscompampl / minmcompampl)**2.
            gdatmodi.nextstdv = stdvstdpcomp / (nextcompampl / minmcompampl)**2.
            gdatmodi.thislfctasym += sum(0.5 * (nextcompunit - thiscompunit)**2 * (1. / gdatmodi.thisstdv**2 - 1. / gdatmodi.nextstdv**2))
        else:
            gdatmodi.thisstdvpara = stdvstdpcomp / (minimum(thiscompampl, nextcompampl) / minmcompampl)**0.5
    else:
        gdatmodi.thisstdvpara = stdvstdpcomp


def show_dbug(gdat, gdatmodi):
    print 'hey'
    print 'thissamp nextsamp diffsampvarb'
    for k in gdat.fittindxpara:
        print '%10.3g  %10.3g %10.3g' % (gdatmodi.thissamp[k], gdatmodi.nextsamp[k], gdatmodi.nextsampvarb[k] - gdatmodi.thissampvarb[k])
    print

        
def retr_indxsamppnts(gdat, l, indxelem):

    indxsamppnts = gdat.fittindxsampcompinit + gdat.fittnumbtrapcuml[l] + indxelem[None, :] * gdat.fittnumbcomp[l] + gdat.fittindxcomp[l][:, None]
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
    psec = (abs(fft.fft2(conv))**2)[:gdat.numbsidecart/2, :gdat.numbsidecart/2] * 1e-3
    
    return psec
   

def retr_psecodim(gdat, psec):

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

                #    plot_genemaps(gdat, gdatmodi, 'this', 'cntpdata', thisindxener=0, thisindxevtt=0, cond=True)
                
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
    

def retr_cntspnts(gdat, lgal, bgal, spec):
    
    indxpixltemp = retr_indxpixl(gdat, bgal, lgal)
    cnts = zeros((gdat.numbener, lgal.size))
    for k in range(lgal.size):
        cnts[:, k] += spec[:, k] * gdat.expototl[:, indxpixltemp[k]]
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

    gdat.lablcurv = r'\kappa'
    gdat.lablexpc = r'E_{c}'
    
    gdat.scalcurvplot = 'self'
    gdat.scalexpcplot = 'self'
    gdat.factcurvplot = 1.
    gdat.factexpcplot = 1.
    
    # temp
    if gdat.elemtype == 'lght':
        gdat.listnamefeateval = ['lgal', 'bgal', 'spec']
        gdat.liststrgfeatplot = []
    if gdat.elemtype == 'lens':
        gdat.listnamefeateval = []
        gdat.liststrgfeatplot = []
    if gdat.elemtype == 'clus':
        gdat.listnamefeateval = ['lgal', 'bgal', 'nobj']
        gdat.liststrgfeatplot = []
    
    # axes
    gdat.minmlgaldata = -gdat.maxmgangdata
    gdat.maxmlgaldata = gdat.maxmgangdata
    gdat.minmbgaldata = -gdat.maxmgangdata
    gdat.maxmbgaldata = gdat.maxmgangdata
    
    # input data
    gdat.pathinpt = gdat.pathdata + 'inpt/'
    if gdat.datatype == 'inpt':
        path = gdat.pathinpt + gdat.strgexprsbrt
        gdat.sbrtdata = pf.getdata(path)
        if gdat.pixltype == 'cart' and not gdat.forccart:
            if gdat.sbrtdata.ndim != 4:
                raise Exception('exprsbrtdata should be a 4D numpy array if pixelization is Cartesian.')
        else:
            if gdat.sbrtdata.ndim != 3:
                raise Exception('exprsbrtdata should be a 3D numpy array if pixelization is HealPix.')
        
        if gdat.pixltype == 'cart':
            if gdat.forccart:
                gdat.sbrtdatatemp = empty((gdat.sbrtdata.shape[0], gdat.numbsidecart, gdat.numbsidecart, gdat.sbrtdata.shape[2]))
                for i in arange(gdat.sbrtdata.shape[0]):
                    for m in arange(gdat.sbrtdata.shape[2]):
                        gdat.sbrtdatatemp[i, :, :, m] = tdpy.util.retr_cart(gdat.sbrtdata[i, :, m], numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).T
                gdat.sbrtdata = gdat.sbrtdatatemp
            else:
                gdat.numbsidecart = gdat.sbrtdata.shape[1]
            gdat.sbrtdata = gdat.sbrtdata.reshape((gdat.sbrtdata.shape[0], gdat.numbsidecart**2, gdat.sbrtdata.shape[3]))
        
        # temp
        #if gdat.strgcnfg.startswith('pcat_ferm'):
        #    gdat.sbrtdata[:, :, 0] = gdat.sbrtdata[:, :, 2]
        #    gdat.sbrtdata[:, :, 1] = gdat.sbrtdata[:, :, 3]
        #    gdat.sbrtdata = gdat.sbrtdata[:, :, 0:2]

        gdat.numbenerfull = gdat.sbrtdata.shape[0]
        gdat.numbpixlfull = gdat.sbrtdata.shape[1]
        gdat.numbevttfull = gdat.sbrtdata.shape[2]
        
        if gdat.pixltype == 'heal':
            # temp
            gdat.numbsidecart = 100
            gdat.numbsideheal = int(sqrt(gdat.numbpixlfull / 12))
        
    gdat.indxenerfull = arange(gdat.numbenerfull)
    gdat.indxevttfull = arange(gdat.numbevttfull)
    
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

    if gdat.elemtype == 'lght':
        gdat.lablelemextn = r'\rm{pts}'
    if gdat.elemtype == 'lens':
        gdat.lablelemextn = r'\rm{sub}'
    if gdat.elemtype == 'clus':
        gdat.lablelemextn = r'\rm{cls}'
    

    # set up the indices of the fitting model
    retr_indxsamp(gdat, init=True)
    
    gdat.listnamefeatrefr = [[] for q in gdat.indxrefr]
    gdat.listnamefeatrefronly = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
    
    # set up the indices of the fitting model
    retr_indxsamp(gdat)
    
    gdat.listnamevarbstat = ['samp', 'sampvarb', 'indxelemfull', 'indxsampcomp', 'lliktotl', 'llik', 'lpritotl']
    for name in gdat.fittlistnamediff:
        #if not (name.startswith('back') and gdat.fittunifback[int(name[4:])]):
        gdat.listnamevarbstat += ['sbrt' + name + 'conv']
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        gdat.listnamevarbstat += ['sbrtpnts']
    if gdat.fittlensmodltype != 'none':
        gdat.listnamevarbstat += ['deflhost']
    if gdat.fittpsfnevaltype == 'conv' or gdat.fittpsfnevaltype == 'full':
        gdat.listnamevarbstat += ['psfnconv']
    if gdat.elemtype == 'lens':
        gdat.listnamevarbstat += ['deflelem']
    
    gdat.liststrgmodl = ['fitt']
    if gdat.datatype == 'mock':
        gdat.liststrgmodl += ['true']
    
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
 
    # temp
    ## names of the variables for which cumulative posteriors will be plotted
    #if gdat.numbswep >= 1e6 and gdat.elemtype == 'lens':
    #    gdat.listnamevarbcpos = ['convelem']
    #else:
    #    gdat.listnamevarbcpos = []
    gdat.listnamevarbcpos = ['convelem']

    if gdat.elemtype == 'lens':
        gdat.ascaglob = 0.05 / gdat.anglfact
        gdat.acutglob = 1. / gdat.anglfact
    
    gdat.cutfdefs = 3e-3 / gdat.anglfact

    if gdat.elemtype == 'lght':
        gdat.namefeateval = 'flux'
    if gdat.elemtype == 'lens':
        gdat.namefeateval = 'defs'
    if gdat.elemtype == 'clus':
        gdat.namefeateval = 'nobj'
    
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
        gdat.lablfracsubh = r'$f_{\rm{sub}}$'

        gdat.lablfracsubhdeltbein = r'$f_{\rm{E,sub}}$'
        gdat.scalfracsubhdeltbein = 'logt'
        
        gdat.lablfracsubhintgbein = r'$f_{\rm{<E,sub}}$'
        gdat.scalfracsubhintgbein = 'logt'
        
        gdat.lablmasssubhdeltbein = r'$M_{\rm{E,sub}}$'
        gdat.scalmasssubhdeltbein = 'logt'
        
        gdat.lablmasssubhintgbein = r'$M_{\rm{<E,sub}}$'
        gdat.scalmasssubhintgbein = 'logt'
        
        gdat.lablmasshostdeltbein = r'$M_{\rm{E,hst}}$'
        gdat.scalmasshostdeltbein = 'logt'
    
        gdat.lablmasshostintgbein = r'$M_{\rm{<E,hst}}$'
        gdat.scalmasshostintgbein = 'logt'
    
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
    gdat.lablmassintg = r'<M>_{<r}'
    gdat.lablmassintgunit = r'$M_{\odot}$'
    gdat.lablmassdelt = r'<M>_r'
    gdat.lablmassdeltunit = r'$M_{\odot}$'
    gdat.lablfracsubhintg = r'<f>_{\rm{<r,sub}}'
    gdat.lablfracsubhdelt = r'<f>_{\rm{r,sub}}'

    gdat.labldefs = r'\alpha_s'
    gdat.lablflux = 'f'
    gdat.lablnobj = 'p'
    
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
    gdat.lablfluxunit = getattr(gdat, 'lablfluxene0' + gdat.nameenerunit + 'unit')
    gdat.lablsbrtunit = getattr(gdat, 'lablsbrtene0' + gdat.nameenerunit + 'sterunit')

    for q in gdat.indxrefr:
        for l in gdat.fittindxpopl:
            setattr(gdat, 'lablcmplref%dpop%d' % (q, l), '$c_{%d%d}$' % (q, l))
            setattr(gdat, 'lablfdisref%dpop%d' % (q, l), '$f_{%d%d}$' % (q, l))
    
    gdat.lablexpo = r'$\epsilon$'
    gdat.lablexpounit = 'cm$^2$ s'
    
    gdat.lablprvl = '$p$'
    
    gdat.lablsind = 's'
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
                if lablunit == '':
                    setattr(gdat, 'labl' + name + nameextn, '%s$^{%s}$ [%s]' % (labl, strgextn, lablunit))
                else:
                    setattr(gdat, 'labl' + name + nameextn, '$%s^{%s}$' % (labl, strgextn))
    
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
    
    for q in gdat.indxrefr:
        setattr(gdat, 'minmcmplpop%d' % q, 0.)
        setattr(gdat, 'maxmcmplpop%d' % q, 1.)
        setattr(gdat, 'corrcmplpop%d' % q, 1.)
        setattr(gdat, 'factcmplpop%dplot' % q, 1.)
        setattr(gdat, 'scalcmplpop%d' % q, 'self')

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

    # feature plotting factors and scalings
    gdat.lablfeat = {}
    gdat.dictglob = {}
    gdat.lablfeatunit = {}
    gdat.lablfeattotl = {}
    for strgmodl in gdat.liststrgmodl:
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        for strgfeat in liststrgfeattotl + gdat.liststrgfeatplot:
            
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
    
    # construct the fitting model
    setp_fixp(gdat)
    
    print 'gdat.fittindxpopl'
    print gdat.fittindxpopl

    for l in gdat.fittindxpopl:
        setattr(gdat, 'minmfdispop%d' % l, 0.)
        setattr(gdat, 'maxmfdispop%d' % l, 1.)
        setattr(gdat, 'corrfdispop%d' % l, 1.)
        setattr(gdat, 'factfdispop%dplot' % l, 1.)
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
    if gdat.fittlensmodltype != 'none' or gdat.truelensmodltype != 'none':
        for strgtemp in ['delt', 'intg']:
            gdat.listnamevarbscal += ['masshost' + strgtemp + 'bein']
            setattr(gdat, 'minmmasshost' + strgtemp + 'bein', gdat.minmmasshost)
            setattr(gdat, 'maxmmasshost' + strgtemp + 'bein', gdat.maxmmasshost)
        if gdat.fittnumbtrap > 0:
            if gdat.elemtype == 'lens':
                for strgtemp in ['delt', 'intg']:
                    gdat.listnamevarbscal += ['masssubh' + strgtemp + 'bein', 'fracsubh' + strgtemp + 'bein'] 
                    setattr(gdat, 'minmmasssubh' + strgtemp + 'bein', gdat.minmmasssubh)
                    setattr(gdat, 'maxmmasssubh' + strgtemp + 'bein', gdat.maxmmasssubh)
                    setattr(gdat, 'minmfracsubh' + strgtemp + 'bein', gdat.minmfracsubh)
                    setattr(gdat, 'maxmfracsubh' + strgtemp + 'bein', gdat.maxmfracsubh)
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
        setattr(gdat, 'fittminmnumbelempop%d' % l, gdat.fittminmnumbelem[l])
        setattr(gdat, 'fittmaxmnumbelempop%d' % l, gdat.fittmaxmnumbelem[l])
    
    # log-prior register
    ## indices of penalization term
    indxlpripena = 0
    ## indices of split and merge term
    indxlprispme = -1
    ## number of elements
    numb = 0
    for l in gdat.fittindxpopl:
        numb += len(gdat.fittliststrgfeat[l])
    if gdat.fittmaxmnumbelemtotl > 0:
        gdat.numblpri = 1 + gdat.fittnumbpopl + numb
    else:
        gdat.numblpri = 0
    if gdat.penalpridiff:
        gdat.numblpri += 1

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
        gdat.mdencrit = retr_mdencrit(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour)
        
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
            # temp
            #gdat.maxmangl = 15. / gdat.anglfact
            gdat.maxmangl = 100. / gdat.anglfact
        if gdat.exprtype == 'chan':
            gdat.maxmangl = 8. / gdat.anglfact
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
        gdat.legdrefr = ['True' for l in gdat.trueindxpopl]
    gdat.legdrefrmiss = []
    gdat.legdrefrhits = []
    for q in gdat.indxrefr:
        gdat.legdrefrmiss.append(gdat.legdrefr[q] + ' miss')
        gdat.legdrefrhits.append(gdat.legdrefr[q] + ' hit')
    if gdat.datatype == 'mock':
        gdat.legdrefrhost = gdat.legdrefr[0] + ' host'
        gdat.legdrefrsour = gdat.legdrefr[0] + ' sour'
   
    # temp
    # element indices in each population
    #gdat.indxelempopl = []
    #for l in gdat.indxpopl:
    #    gdat.indxelempopl.append(arange(sum(gdat.maxmnumbelem[:l]), sum(gdat.maxmnumbelem[:l+1])))
    
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
        gdat.listnamechro += ['pixlkern', 'kernelem']
        gdat.listlegdchro += ['Pixels Kern. Ev.', 'Kern. Ev.']
    if gdat.elemtype == 'lens':
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

    if gdat.elemtype == 'lens': 
        retr_axis(gdat, 'defl', -gdat.maxmgang, gdat.maxmgang, 40)
        retr_axis(gdat, 'deflelem', -gdat.maxmgang * 1e-2, gdat.maxmgang * 1e-2, 40)

    # lensing problem setup
    ## number of deflection components to plot

    gdat.binslgalcart = linspace(gdat.minmlgaldata, gdat.maxmlgaldata, gdat.numbsidecart + 1)
    gdat.binsbgalcart = linspace(gdat.minmbgaldata, gdat.maxmbgaldata, gdat.numbsidecart + 1)
    gdat.meanlgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
    gdat.meanbgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
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
        gdat.meanmpollgal, gdat.meanmpolbgal = meshgrid(gdat.meanmpolodim, gdat.meanmpolodim, indexing='ij')
        gdat.meanmpol = sqrt(gdat.meanmpollgal**2 + gdat.meanmpolbgal**2)

    # element parameter vector indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompflux = 2
    gdat.indxcompsind = 3
    gdat.indxcompcurv = 4
    gdat.indxcompexpc = 4
    # temp
    #gdat.indxcomp = [[] for l in gdat.indxpopl]
    #for l in gdat.indxpopl:
    #    gdat.indxcomp[l] = arange(gdat.fittnumbcomp[l])
    #gdat.indxelem = []
    #for l in gdat.indxpopl:
    #    gdat.indxelem.append(arange(gdat.maxmnumbelem[l]))

    # exposure
    if gdat.correxpo:
        if isinstance(gdat.strgexpo, float):
            if gdat.datatype == 'mock':
                if gdat.pixltype == 'heal':
                    gdat.expo = gdat.strgexpo * ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                if gdat.pixltype == 'cart':
                    gdat.expo = gdat.strgexpo * ones((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull))
            if gdat.datatype == 'inpt':
                gdat.expo = gdat.strgexpo * ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
        else: 
            path = gdat.pathinpt + gdat.strgexpo
            gdat.expo = pf.getdata(path)
            if amin(gdat.expo) == amax(gdat.expo):
                raise Exception('Bad input exposure map.')
            if gdat.pixltype == 'cart':
   
                if gdat.forccart:
                    gdat.expotemp = empty((gdat.numbenerfull, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevttfull))
                    for i in gdat.indxenerfull:
                        for m in gdat.indxevttfull:
                            gdat.expotemp[i, :, :, m] = tdpy.util.retr_cart(gdat.expo[i, :, m], numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                                                             minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                                             minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).T
                    gdat.expo = gdat.expotemp
                
                gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))

    if gdat.killexpo:
        gdat.expo *= 1e-90

    # only include desired energy and PSF class bins 
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
    ## exposure
    if gdat.correxpo:
        gdat.expo = gdat.expo[gdat.indxcubeincl]
    ## backgrounds
    for strgmodl in gdat.liststrgmodl:
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        indxback = getattr(gdat, strgmodl + 'indxback')
        indxtessincl = meshgrid(indxback, gdat.indxenerincl, gdat.indxpixlfull, gdat.indxevttincl, indexing='ij')
        sbrtbacknorm = sbrtbacknorm[indxtessincl]
        setattr(gdat, strgmodl + 'sbrtbacknorm', sbrtbacknorm)
    
    if gdat.datatype == 'inpt':
        gdat.sbrtdata = gdat.sbrtdata[gdat.indxcubeincl]
        gdat.cntpdata = retr_cntp(gdat, gdat.sbrtdata)

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
    if gdat.mask != None:
        if gdat.masktype == 'sqre':
            indxpixlmask = where((gdat.lgalgrid > gdat.mask[0]) & (gdat.lgalgrid < gdat.mask[1]) & (gdat.bgalgrid > gdat.mask[2]) & (gdat.bgalgrid < gdat.mask[3]))[0]
            gdat.expo[:, indxpixlmask, :] = 0.

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
    
    # temp
    gdat.expo[where(gdat.expo < 1e-50)] = 1e-50
    
    # exclude voxels with vanishing exposure
    if gdat.correxpo:
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                gdat.indxpixlrofi = intersect1d(gdat.indxpixlrofi, where(gdat.expo[i, :, m] > 0.)[0])
   
    gdat.minmexpo = amin(gdat.expo[where(gdat.expo > 1e-100)])
    gdat.maxmexpo = amax(gdat.expo)

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
    
    gdat.indxcuberofi = meshgrid(gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
    gdat.numbpixl = gdat.indxpixlrofi.size
    gdat.indxpixl = arange(gdat.numbpixl)
    gdat.indxcube = meshgrid(gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
   
    ## meshed indcices for the four-dimensional background data structures
    for strgmodl in gdat.liststrgmodl:
        indxback = getattr(gdat, strgmodl + 'indxback')
        indxtessback = []
        for c in indxback:  
            indxtessback.append(meshgrid(array([c]), gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij'))
        setattr(gdat, strgmodl + 'indxtessback', indxtessback)
        indxtess = meshgrid(indxback, gdat.indxener, gdat.indxpixl, gdat.indxevtt, indexing='ij')
        setattr(gdat, strgmodl + 'indxtess', indxtess)

    ## exposure
    if gdat.correxpo:
        gdat.expofull = copy(gdat.expo)
        gdat.expo = gdat.expo[gdat.indxcuberofi]
    ## backgrounds
    for strgmodl in gdat.liststrgmodl:
        indxback = getattr(gdat, strgmodl + 'indxback')
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        
        if gdat.pixltype == 'heal':
            sbrtbackhealfull = copy(sbrtbacknorm[meshgrid(indxback, gdat.indxener, gdat.indxpixlfull, gdat.indxevtt, indexing='ij')])
            setattr(gdat, strgmodl + 'sbrtbackhealfull', sbrtbackhealfull)
        indxtessrofi = meshgrid(indxback, gdat.indxener, gdat.indxpixlrofi, gdat.indxevtt, indexing='ij')
        sbrtbacknorm = sbrtbacknorm[indxtessrofi]
        setattr(gdat, strgmodl + 'sbrtbacknorm', sbrtbacknorm)
                
    gdat.expototl = sum(gdat.expo, axis=2)
    gdat.expototlmean = mean(gdat.expototl, axis=1)

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
            gdat.pixlcnvt = zeros(gdat.numbpixl, dtype=int) - 1
            numbpixlmarg = gdat.indxpixlrofimarg.size
            for k in range(numbpixlmarg):
                dist = retr_angldistunit(gdat, lgalheal[gdat.indxpixlrofimarg[k]], bgalheal[gdat.indxpixlrofimarg[k]], gdat.indxpixl)
                gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
    
    # dummy pixel indices for full (nonlocal) element kernel evaluation 
    gdat.listindxpixl = []
    numb = max(sum(gdat.fittmaxmnumbelem), sum(gdat.truemaxmnumbelem)) + 2
    for k in range(numb):
        gdat.listindxpixl.append([])
        for kk in range(k):
            gdat.listindxpixl[k].append(gdat.indxpixl)

    if gdat.evalcirc != 'full' and gdat.numbpixl * gdat.fittmaxmnumbelemtotl < 1e5:
        # temp
        #gdat.calcerrr = True
        gdat.calcerrr = False
    else:
        gdat.calcerrr = False
   
    # spatial averaging setup
    gdat.listnamespatmean = ['full']
    gdat.listindxcubespatmean = [gdat.indxcube]
    gdat.numbspatmean = len(gdat.listnamespatmean)

    gdat.psfnexpr = retr_psfn(gdat, gdat.psfpexpr, gdat.indxener, gdat.binsangl, gdat.psfntypeexpr, gdat.binsoaxi, gdat.exproaxitype)
    
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
                    psfnwdth = retr_psfnwdth(gdat, gdat.psfnexpr, frac)
                    gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                    gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
        
        if gdat.exprtype == 'chan':
            gdat.maxmangleval = maximum(gdat.maxmangleval, array([3., 5., 7.]) / gdat.anglfact)
       
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
                        if gdat.maxmangl < sqrt(2.) * gdat.maxmgang:
                            raise Exception('Angular axis used to interpolate the PSF should be longer.')
                    gdat.indxpixlprox[h].append(indxpixlproxtemp)
                cntrsave = tdpy.util.show_prog(j, gdat.numbpixl, cntrsave)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
        gdat.numbpixlprox = zeros(gdat.numbprox) 
        for k in range(gdat.numbprox):
            for m in range(len(gdat.indxpixlprox[k])):
                if isinstance(gdat.indxpixlprox[k][m], int):
                    gdat.numbpixlprox[k] += gdat.numbpixl
                else:
                    gdat.numbpixlprox[k] += len(gdat.indxpixlprox[k][m])
            gdat.numbpixlprox[k] /= len(gdat.indxpixlprox[k])
        print 'gdat.numbpixlprox'
        print gdat.numbpixlprox
    else:
        gdat.numbpixlprox = gdat.numbpixl

    # turn off relevant proposal types
    gdat.indxfixpprop = []
    for k, strg in enumerate(gdat.fittnamefixp):
        if gdat.fittnumbtrap > 0 and k in gdat.fittindxfixpnumbelem:
            thisbool = False
        else:
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
    gdat.propfixp = gdat.propmeanelem or gdat.propbacp or gdat.proppsfp or gdat.proplenp
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
            gdat.indxstdpexpc = gdat.numbfixpprop + 4
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
    indxelemfull = [range(gdat.fittmaxmnumbelem[l]) for l in gdat.fittindxpopl]
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
        gdat.stdvstdp = 1e-4 + zeros(gdat.numbstdp)

        if gdat.fittnumbtrap > 0:
            gdat.stdvstdp[gdat.indxstdpmeanelempop0] = 1e-1
            gdat.stdvstdp[gdat.indxstdpdefsdistsloppop0] = 1e-1

        gdat.stdvstdp[gdat.indxstdpsigcene0evt0] = 3e-2
        gdat.stdvstdp[gdat.indxstdpbacpbac0ene0] = 1e-3
        
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
                indxback = where(gdat.fittnameback == 'isot')[0]
                if indxback.size > 0:
                    gdat.stdvstdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpbac%dene%d' % (indxback[0], i))]] = 5e-3
            for i in gdat.indxener:
                indxback = where(gdat.fittnameback == 'fdfm')[0]
                if indxback.size > 0:
                    gdat.stdvstdp[gdat.indxstdppara[getattr(gdat, 'fittindxfixpbacpbac%dene%d' % (indxback[0], i))]] = 5e-4
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

    ## input data
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
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
    gdat.numbpixlsave = min(10000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)
    gdat.indxcubesave = meshgrid(gdat.indxener, gdat.indxpixlsave, gdat.indxevtt, indexing='ij')
   
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
            
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        gdat.exprfwhm = 2. * retr_psfnwdth(gdat, gdat.psfnexpr, 0.5)
        gdat.stdvspatprio = amax(gdat.exprfwhm)
    if gdat.elemtype == 'lens':
        gdat.stdvspatprio = amax(gdat.psfpexpr)
    
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
    
    ### filters
    #### filter for reference elements
    gdat.listnametruefilt = ['']
    gdat.listindxtrueelemfilt = []
    for nametruefilt in gdat.listnametruefilt:
        indxtrueelemfilt = [[] for q in gdat.indxrefr]
        if nametruefilt == '':
            for q in gdat.indxrefr:
                indxtrueelemfilt[q] = arange(gdat.truenumbelem[q])
        gdat.listindxtrueelemfilt.append(indxtrueelemfilt)

    #### filter for model elements
    gdat.listnamemodlfilt = ['']
    if gdat.priofactdoff <= 1.:
        gdat.listnamemodlfilt += ['deltllik']
    #### model elements inside the image
    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
        gdat.listnamemodlfilt += ['imag']
    #### model subhalos inside high normalized relevance region
    if gdat.elemtype == 'lens':
        gdat.listnamemodlfilt += ['nrel']
    
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
                    
    if gdat.datatype == 'inpt':
        proc_cntpdata(gdat)
    
    gdat.sbrtpntstemp = empty_like(gdat.expo)


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
                            temp = retr_rele(gdat, gdat.truecntpdiffconv['lens'][0, :, 0], gdat.binslgal[k], gdat.binsbgal[n], 1., \
                                                                                                    gdat.trueascadistmeanpop0, gdat.trueacutdistmeanpop0, gdat.indxpixl)
    
                            #temp /= amax(temp)
                            pdfnpriotemp[k, n] = temp**gdat.relnpowr
                    lpdfprio, lpdfprioobjt = retr_spatprio(gdat, pdfnpriotemp)
                    lpdfpriointp = lpdfprioobjt(gdat.meanbgalcart, gdat.meanlgalcart)
        setattr(gdat, strgmodl + 'lpdfprio', lpdfprio)
        setattr(gdat, strgmodl + 'lpdfprioobjt', lpdfprioobjt)
        setattr(gdat, strgmodl + 'lpdfpriointp', lpdfpriointp)
        
    # data structures for JIT
    gdat.deflelem = zeros((gdat.numbpixl, 2))

    # plot settings
    ## upper limit of histograms
    if gdat.datatype == 'inpt':
        gdat.limtpntshist = [0.5, 10**ceil(log10(gdat.fittmaxmnumbelemtotl))]
   
    ## spatial average
    gdat.sbrtdatamean = retr_spatmean(gdat, gdat.cntpdata, boolcntp=True)
    
    # find the pixels at which data count maps have local maxima
    if gdat.pixltype == 'cart':
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                # temp
                gdat.indxxdatmaxm, gdat.indxydatmaxm = tdpy.util.retr_indximagmaxm(gdat.cntpdatacart[i, :, m])

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
    grad = dstack((gradient(maps, gdat.sizepixl, axis=0), gradient(maps, gdat.sizepixl, axis=1)))

    return grad


def retr_spatmean(gdat, inpt, boolcntp=False):
    
    listspatmean = empty((gdat.numbspatmean, gdat.numbener))
    for b, namespatmean in enumerate(gdat.listnamespatmean):
        if boolcntp:
            cntp = inpt[gdat.listindxcubespatmean[b]]
        else:
            if False and gdat.exprtype == 'hubb':
                print 'inpt'
                summgene(inpt)
                print 'gdat.expo'
                summgene(gdat.expo)
                print 'gdat.listindxcubespatmean[0]'
                summgene(gdat.listindxcubespatmean[b][0])
                print 'gdat.listindxcubespatmean[1]'
                summgene(gdat.listindxcubespatmean[b][1])
                print 'gdat.listindxcubespatmean[2]'
                summgene(gdat.listindxcubespatmean[b][2])
            cntp = inpt[gdat.listindxcubespatmean[b]] * gdat.expo[gdat.listindxcubespatmean[b]] * gdat.apix
            if gdat.enerdiff:
                cntp *= gdat.deltener[:, None, None]
        
        spatmean = mean(sum(cntp, 2), axis=1) / gdat.expototlmean / gdat.apix
        if gdat.enerdiff:
            spatmean /= gdat.deltener
        
        listspatmean[b, :] = spatmean

    return listspatmean


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
    gdat.maxmcntpdata = max(0.5, amax(sum(gdat.cntpdata, 2)))
    gdat.minmcntpdata = max(0.5, 1e-3 * gdat.maxmcntpdata)
    gdat.maxmcntpresi = ceil(gdat.maxmcntpdata * 0.1)
    gdat.minmcntpresi = -gdat.maxmcntpresi

    # plotting
    liststrgcbar = ['cntpdata', 'cntpresi']
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
    

def retr_fromgdat(gdat, gdatmodi, strg, strgvarb, mometype='medi', indxvarb=None, indxlist=None):
    
    if strgvarb == 'cntpdata':
        varb = gdat.cntpdata
    else:
        if strg == 'true':
            varb = getattr(gdat, strg + strgvarb)
        if strg == 'this':
            if mometype == 'errr':
                varb = getattr(gdatmodi, strg + 'errr' + strgvarb)
            else:
                varb = getattr(gdatmodi, strg + strgvarb)
        if strg == 'post':
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
    if init:
        # transdimensional element populations
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        indxpopl = arange(numbpopl, dtype=int) 
        
        numbback = len(getattr(gdat, strgmodl + 'backtype'))
        indxback = arange(numbback)
    
    else:
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        indxback = getattr(gdat, strgmodl + 'indxback')
        hostemistype = getattr(gdat, strgmodl + 'hostemistype')
        maxmnumbelem = getattr(gdat, strgmodl + 'maxmnumbelem') 
        
        # construct background surface brightness templates from the user input
        lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
        backtype = getattr(gdat, strgmodl + 'backtype')
        numbback = getattr(gdat, strgmodl + 'numbback')
        sbrtbacknorm = empty((numbback, gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
        unifback = ones(numbback, dtype=bool)
        for c in indxback:
            if backtype[c] == 'data':
                sbrtbacknormtemp = copy(gdat.sbrtdata)
                sbrtbacknormtemp[where(sbrtbacknormtemp == 0.)] = 1e-100
            elif isinstance(backtype[c], float):
                if gdat.pixltype == 'heal':
                    sbrtbacknormtemp = zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + backtype[c]
                if gdat.pixltype == 'cart':
                    sbrtbacknormtemp = zeros((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + backtype[c]
            elif isinstance(backtype[c], ndarray) and backtype[c].ndim == 1:
                if gdat.pixltype == 'heal':
                    sbrtbacknormtemp = zeros((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull)) + backtype[c][:, None, None]
                if gdat.pixltype == 'cart':
                    sbrtbacknormtemp = zeros((gdat.numbenerfull, gdat.numbsidecart**2, gdat.numbevttfull)) + backtype[c][:, None, None]
            else:
                path = gdat.pathinpt + backtype[c]
                sbrtbacknormtemp = pf.getdata(path)
                
                if gdat.pixltype == 'cart':
                    if gdat.forccart:
                        sbrtbacknormtemptemp = empty((gdat.numbenerfull, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevttfull))
                        for i in gdat.indxenerfull:
                            for m in gdat.indxevttfull:
                                sbrtbacknormtemptemp[i, :, :, m] = tdpy.util.retr_cart(sbrtbacknormtemp[i, :, m], numbsidelgal=gdat.numbsidecart, numbsidebgal=gdat.numbsidecart, \
                                                                                                 minmlgal=gdat.anglfact*gdat.minmlgaldata, maxmlgal=gdat.anglfact*gdat.maxmlgaldata, \
                                                                                                 minmbgal=gdat.anglfact*gdat.minmbgaldata, maxmbgal=gdat.anglfact*gdat.maxmbgaldata).T
                        sbrtbacknormtemp = sbrtbacknormtemptemp
                    
                    if sbrtbacknormtemp.shape[1] != gdat.numbsidecart:
                        print 'gdat.numbsidecart'
                        print gdat.numbsidecart
                        print 'sbrtbacknormtemp.shape[1]'
                        print sbrtbacknormtemp.shape[1]
                        raise Exception('Provided background template must have the chosen image dimensions.')
            
                    sbrtbacknormtemp = sbrtbacknormtemp.reshape((sbrtbacknormtemp.shape[0], -1, sbrtbacknormtemp.shape[-1]))

            # temp
            if gdat.strgcnfg.startswith('pcat_ferm') and sbrtbacknormtemp.shape[-1] == 2:
                sbrtbacknorm[c, :, :, :] = 1.
                sbrtbacknorm[c, :, :, 2:4] = sbrtbacknormtemp
            else:
                sbrtbacknorm[c, ...] = sbrtbacknormtemp
            
            # determine spatially uniform background templates
            for i in gdat.indxenerfull:
                for m in gdat.indxevttfull:
                    if std(sbrtbacknorm[c, i, :, m]) > 1e-6:
                        unifback[c] = False

            if amin(sbrtbacknorm[c, ...]) < 0.:
                booltemp = False
                raise Exception('Background templates must be positive-definite everywhere.')
        
        boolzero = True
        for c in indxback:
            if amin(sbrtbacknorm[c, ...]) > 0. or backtype[c] == 'data':
                boolzero = False
        if boolzero:
            raise Exception('At least one background template must be positive everywhere.')
       
        if 'data' in backtype and len(backtype) != 1:
            raise Exception('data - PS residual can be the only background.')
        if backtype[0] != 'data' and gdat.penalpridiff:
            raise Exception('Diffuse background power spectrum penalization is unncessary if the background is not the data - PS residual.')
        
        if sum(maxmnumbelem) > 0 and (gdat.elemtype == 'lght' or gdat.elemtype == 'clus'):
            if hostemistype != 'none' or not unifback.all():
                psfnevaltype = 'full'
            else:
                psfnevaltype = 'kern'
        else:
            if hostemistype != 'none' or not unifback.all():
                psfnevaltype = 'conv'
            else:
                psfnevaltype = 'none'
        
        setp_namevarbvalu(gdat, 'psfnevaltype', psfnevaltype, strgmodl=strgmodl)
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
                        setp_namevarblimt(gdat, 'sigcene%devt%d' % (i, m), [meansigc, stdvsigc], typelimt='meanstdv', strgmodl=strgmodl)
                        if psfntype == 'doubking' or psfntype == 'singking':
                            meangamc = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 1]
                            stdvgamc = meangamc * 0.1
                            setp_namevarblimt(gdat, 'gamcene%devt%d' % (i, m), [meangamc, stdvgamc], typelimt='meanstdv', strgmodl=strgmodl)
                            if psfntype == 'doubking':
                                meansigt = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 2]
                                stdvsigt = meansigt * 0.1
                                setp_namevarblimt(gdat, 'sigtene%devt%d' % (i, m), [meansigt, stdvsigt], typelimt='meanstdv', strgmodl=strgmodl)
                                meangamt = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 3]
                                stdvgamt = meangamt * 0.1
                                setp_namevarblimt(gdat, 'gamtene%devt%d' % (i, m), [meangamt, stdvgamt], typelimt='meanstdv', strgmodl=strgmodl)
                                meanpsff = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 4]
                                stdvpsff = meanpsff * 0.1
                                setp_namevarblimt(gdat, 'psffene%devt%d' % (i, m), [meanpsff, stdvpsff], typelimt='meanstdv', strgmodl=strgmodl)
                        elif oaxitype:
                            meanonor = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 1]
                            stdvonor = meanonor * 0.1
                            setp_namevarblimt(gdat, 'onorene%devt%d' % (i, m), [meanonor, stdvonor], typelimt='meanstdv', strgmodl=strgmodl)
                            meanoind = gdat.psfpexpr[i * numbpsfpform + m * numbpsfpform * gdat.numbener + 2]
                            stdvoind = meanoind * 0.1
                            setp_namevarblimt(gdat, 'oindene%devt%d' % (i, m), [meanoind, stdvoind], typelimt='meanstdv', strgmodl=strgmodl)
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
                setp_namevarblimt(gdat, 'sigc', [minmsigm, maxmsigm], ener=True, evtt=True, strgmodl=strgmodl)
                setp_namevarblimt(gdat, 'sigt', [minmsigm, maxmsigm], ener=True, evtt=True, strgmodl=strgmodl)
                setp_namevarblimt(gdat, 'gamc', [minmgamm, maxmgamm], ener=True, evtt=True, strgmodl=strgmodl)
                setp_namevarblimt(gdat, 'gamt', [minmgamm, maxmgamm], ener=True, evtt=True, strgmodl=strgmodl)
                setp_namevarblimt(gdat, 'onor', [minmonor, maxmonor], ener=True, evtt=True, strgmodl=strgmodl)
                setp_namevarblimt(gdat, 'oind', [minmoind, maxmoind], ener=True, evtt=True, strgmodl=strgmodl)
            setp_namevarblimt(gdat, 'psff', [0., 1.], ener=True, evtt=True, strgmodl=strgmodl)
 
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
        indxbacpback = []
        indxbackbacp = zeros(numbbacp, dtype=int)
        indxenerbacp = zeros(numbbacp, dtype=int)
        cntr = 0
        for c in indxback: 
            
            if specback[c]:
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
        maxmnumbelemtotl = sum(maxmnumbelem)
        indxelemtotl = arange(maxmnumbelemtotl)
        maxmnumbelemcumr = cumsum(maxmnumbelem)
        maxmnumbelemcuml = concatenate((array([0]), maxmnumbelemcumr[:-1]))
   
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
                if spectype[l] == 'expc':
                    liststrgcomp[l] += ['expc']
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
                if spectype[l] == 'expc':
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
            if gdat.elemtype == 'lght':
                liststrgfeatodim[l] += ['cnts']
            if gdat.elemtype == 'lens':
                liststrgfeatodim[l] += ['mcut', 'diss', 'rele', 'reln', 'reld', 'relc']
            # temp
            if strgmodl == 'true':
                for namefeat in gdat.listnamefeatsele:
                    for namesele in gdat.listnamesele:
                        liststrgfeatodim[l] += [namefeat + namesele]
        
        # add reference element features that are not available in the PCAT element model
        if strgmodl == 'fitt':
            # temp
            for name, varb in gdat.__dict__.iteritems():
                if name.startswith('refr') and name != 'refrinfo':
                    for q in gdat.indxrefr: 
                        nametemp = name[4:] + gdat.listnamerefr[q]
                        if not nametemp in gdat.listnamefeatrefr[q]:
                            gdat.listnamefeatrefr[q].append(name[4:])
                        for l in indxpopl:
                            if not name[4:] in liststrgfeatodim[l] and name[4:] != 'spec' and name[4:] != 'deflprof' and name[4:] != 'specplot':
                                liststrgfeatodim[l].append(name[4:])
                                if not name[4:] in gdat.listnamefeatrefronly[q][l]:
                                    gdat.listnamefeatrefronly[q][l].append(name[4:])

        if strgmodl == 'true':
            gdat.listnamefeatrefr = liststrgfeatodim

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
        
        # temp
        for l in indxpopl:
            if gdat.elemtype == 'lght':
                liststrgfeat[l] += ['spec', 'specplot']
            if gdat.elemtype == 'lens':
                liststrgfeat[l] += ['deflprof']
        
        # variables for which pair-correlations will be plotted
        liststrgfeatcorr = [[] for l in indxpopl]
        if gdat.plotelemcorr:
            for l in indxpopl:
                for strgfeat in liststrgfeatodim[l] + gdat.liststrgfeatplot:
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
        numbtrappopl = maxmnumbelem * numbcomp
        numbtrapcumr = cumsum(numbtrappopl)
        numbtrapcuml = concatenate((array([0]), numbtrapcumr[:-1]))
        numbtrap = sum(numbtrappopl)
        
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
        if gdat.elemtype == 'lght':
            for strgfeat in ['sind', 'curv', 'expc']:
                if not strgfeat in liststrgfeatdefa:
                    liststrgfeatdefa.append(strgfeat)

        if numbtrap > 0:

            # number of elements
            for l in indxpopl:
                dicttemp['indxfixpnumbelempop%d' % l] = cntr.incr()
            
            # hyperparameters
            ## mean number of elements
            for l in indxpopl:
                dicttemp['indxfixpmeanelempop%d' % l] = cntr.incr()
            
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
        
            dicttemp['indxfixpnumbelem'] = []
            dicttemp['indxfixpmeanelem'] = []
            for strg, valu in dicttemp.iteritems():
                if strg[8:].startswith('numbelemp'):
                    dicttemp['indxfixpnumbelem'].append(valu)
                if strg[8:].startswith('meanelemp'):
                    dicttemp['indxfixpmeanelem'].append(valu)
            dicttemp['indxfixpnumbelem'] = array(dicttemp['indxfixpnumbelem'])
            dicttemp['indxfixpmeanelem'] = array(dicttemp['indxfixpmeanelem'])
        
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
                indx = cntr.incr()
                dicttemp['indxfixpbacpbac%d' % c] = indx
                dicttemp['indxfixpbacp'].append(indx)
            else:
                for i in gdat.indxener:
                    indx = cntr.incr()
                    dicttemp['indxfixpbacpbac%dene%d' % (c, i)] = indx
                    dicttemp['indxfixpbacp'].append(indx)
                    
        if False:
            for attr, valu in dicttemp.iteritems():
                if attr.startswith('indxfixpbac'):
                    print attr
                    print valu
                    print

        dicttemp['indxfixpbacp'] = array(dicttemp['indxfixpbacp'])
        dicttemp['indxfixphost'] = []
        dicttemp['indxfixpsour'] = []
        dicttemp['indxfixplenp'] = []
        
        # temp
        #dicttemp['indxfixpanglsour'] = []
        #dicttemp['indxfixpanglhost'] = []
        #dicttemp['indxfixpangllens'] = []
        
        dicttemp['indxfixpspecsour'] = []
        dicttemp['indxfixpspechost'] = []

        if gdat.elemtype == 'lens':
            if lensmodltype != 'none':
                dicttemp['indxfixplgalsour'] = cntr.incr()
                dicttemp['indxfixpbgalsour'] = cntr.incr()
                dicttemp['indxfixpfluxsour'] = cntr.incr()
                if gdat.numbener > 1:
                    dicttemp['indxfixpsindsour'] = cntr.incr()
                dicttemp['indxfixpsizesour'] = cntr.incr()
                dicttemp['indxfixpellpsour'] = cntr.incr()
                dicttemp['indxfixpanglsour'] = cntr.incr()
            if hostemistype != 'none':
                dicttemp['indxfixplgalhost'] = cntr.incr()
                dicttemp['indxfixpbgalhost'] = cntr.incr()
                dicttemp['indxfixpfluxhost'] = cntr.incr()
                if gdat.numbener > 1:
                    dicttemp['indxfixpsindhost'] = cntr.incr()
                dicttemp['indxfixpsizehost'] = cntr.incr()
            if lensmodltype != 'none':
                dicttemp['indxfixpbeinhost'] = cntr.incr()
            if hostemistype != 'none':
                dicttemp['indxfixpellphost'] = cntr.incr()
                dicttemp['indxfixpanglhost'] = cntr.incr()
                dicttemp['indxfixpserihost'] = cntr.incr()
            if lensmodltype != 'none':
                dicttemp['indxfixpsherextr'] = cntr.incr()
                dicttemp['indxfixpsangextr'] = cntr.incr()
                dicttemp['indxfixpsour'] = []

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
                    
                if strg[8:] == 'fluxsour' or strg[8:] == 'sindsour':
                    dicttemp['indxfixpspecsour'].append(valu)

                if strg[8:] == 'fluxhost' or strg[8:] == 'sindhost':
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
                if spectype[l] == 'expc':
                    liststrgspep[l] += ['expc']
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
    
    if gdat.verbtype > 0:
        if strgmodl == 'true':
            strgtemp = 'true'
        if strgmodl == 'fitt':
            strgtemp = 'fitting'
        print 'Building the %s model...' % strgtemp
    listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
    liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
    
    psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')
    
    numbcomp = getattr(gdat, strgmodl + 'numbcomp')
    numbtrapcuml = getattr(gdat, strgmodl + 'numbtrapcuml')
    numbtrapcumr = getattr(gdat, strgmodl + 'numbtrapcumr')
    indxfixp = getattr(gdat, strgmodl + 'indxfixp')
    numbfixp = getattr(gdat, strgmodl + 'numbfixp')
    numbback = getattr(gdat, strgmodl + 'numbback')
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    numbpara = getattr(gdat, strgmodl + 'numbpara')
    nameback = getattr(gdat, strgmodl + 'nameback')
    legdback = getattr(gdat, strgmodl + 'legdback')

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
        if strg[:-1].endswith('pop'):
            
            if numbpopl == 1:
                strgpopl = ''
                strgpoplcomm = ''
            else:
                strgpopl = '%s' % strg[-1]
                strgpoplcomm = ',%s' % strg[-1]

            if namefixp[k].startswith('numbelem'):
                lablfixp[k] = '$N_{%s%s}$' % (gdat.lablelemextn, strgpoplcomm)
                scalfixp[k] = 'pois'
                
            if namefixp[k].startswith('meanelem'):
                lablfixp[k] = r'$\mu_{%s%s}$' % (gdat.lablelemextn, strgpoplcomm)
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
            
            # index of the background parameter
            c = int(strgvarb[7])#indxbackbacp[k-indxfixpbacpinit]

            name = 'bacpbac%d' % c
            
            if specback[c]:
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
            
            #namefixp[k] = name
            lablfixp[k] = '$A_{%s%s}$' % (strgbacktemp, strgenertemp)
            lablfixpunit[k] = gdat.lablsbrtunit
            
            try:
                scalfixp[k] = getattr(gdat, strgmodl + 'scal' + name)
                if gdat.verbtype > 0:
                    print 'Received custom scaling for %s: %s' % (name, scalfixp[k])
            except:
                scalfixp[k] = 'logt'
        
        if gdat.elemtype == 'lens':
            if k in getattr(gdat, strgmodl + 'indxfixplenp'):
                if strgvarb == 'lgalsour':
                    lablfixp[k] = r'$\theta_{\rm{1,src}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'bgalsour':
                    lablfixp[k] = r'$\theta_{\rm{2,src}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'fluxsour':
                    lablfixp[k] = r'$f_{\rm{src}}$'
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if gdat.numbener > 1:
                    if strgvarb == 'sindsour':
                        lablfixp[k] = r'$s_{\rm{src}}$'
                        scalfixp[k] = 'self'
                if strgvarb == 'sizesour':
                    lablfixp[k] = r'$\theta_{\rm{e,src}}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'ellpsour':
                    lablfixp[k] = r'$\epsilon_{\rm{src}}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'anglsour':
                    lablfixp[k] = r'$\phi_{\rm{src}}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'lgalhost':
                    lablfixp[k] = r'$\theta_{1,\rm{hst}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'bgalhost':
                    lablfixp[k] = r'$\theta_{2,\rm{hst}}$'
                    scalfixp[k] = 'self'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'fluxhost':
                    lablfixp[k] = r'$f_{\rm{hst}}$'
                    scalfixp[k] = 'logt'
                    lablfixpunit[k] = gdat.lablfluxunit
                if gdat.numbener > 1:
                    if strgvarb == 'sindhost':
                        lablfixp[k] = r'$s_{\rm{hst}}$'
                        scalfixp[k] = 'self'
                if strgvarb == 'sizehost':
                    lablfixp[k] = r'$\theta_{\rm{e,hst}}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'beinhost':
                    lablfixp[k] = r'$\theta_{\rm{E,hst}}$'
                    scalfixp[k] = 'logt'
                    factfixpplot[k] = gdat.anglfact
                    lablfixpunit[k] = gdat.lablgangunit
                if strgvarb == 'ellphost':
                    lablfixp[k] = r'$\epsilon_{\rm{hst}}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'anglhost':
                    lablfixp[k] = r'$\phi_{\rm{hst}}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'serihost':
                    lablfixp[k] = r'$n_{\rm{S,hst}}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'sherextr':
                    lablfixp[k] = r'$\gamma_{ext}$'
                    scalfixp[k] = 'self'
                if strgvarb == 'sangextr':
                    lablfixp[k] = r'$\phi_{ext}$'
                    scalfixp[k] = 'self'
        
        if scalfixp[k] == 'pois' or scalfixp[k] == 'self' or scalfixp[k] == 'logt' or scalfixp[k] == 'atan':
            if strgvarb.startswith('numbelem'):
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
    listlegdsbrt = deepcopy(legdback)
    if gdat.elemtype == 'lght':
        listlegdsbrt.append('PS')
    if gdat.elemtype == 'lens':
        listlegdsbrt.append('Source')
        listlegdsbrt.append('Host')
    if gdat.elemtype == 'clus':
        listlegdsbrt.append('Uniform')
    listlegdsbrtspec = deepcopy(listlegdsbrt)
    listlegdsbrtspec += ['Data']
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
        indxpopltemp = argmin(where(k // numbtrapcumr == 0))
        indxcomptemp = (k - numbtrapcuml[indxpopltemp]) % numbcomp[indxpopltemp]
        namepara[numbfixp+k] = liststrgcomp[indxpopltemp][indxcomptemp]
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
            setattr(gdat, strgmodl + strg, valu)


def setp_namevarbsing(gdat, strgvarb, popl, ener, evtt, back):
    
    liststrgvarb = []
    if popl:
        for q in gdat.trueindxpopl:
            liststrgvarb.append(strgvarb + 'pop%d' % q)
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
                meanpara = listvalu[0]
                setattr(gdat, strgmodl + 'mean' + strgvarb, meanpara)
            
            try:
                stdvpara = getattr(gdat, strgmodl + 'stdv' + strgvarb)
                if stdvpara == None:
                    raise
            except:
                stdvpara = listvalu[1]
                setattr(gdat, strgmodl + 'stdv' + strgvarb, stdvpara)

            # set minimum and maximum for Gaussian distributed variables
            try:
                getattr(gdat, strgmodl + 'minm' + strgvarb)
            except:
                setattr(gdat, strgmodl + 'minm' + strgvarb, meanpara - stdvpara * 5)
 
            try:
                getattr(gdat, strgmodl + 'maxm' + strgvarb)
            except:
                setattr(gdat, strgmodl + 'maxm' + strgvarb, meanpara + stdvpara * 5)
 

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


def init_figr(gdat, gdatmodi, strgplot, strg, indxenerplot=None, indxevttplot=None, indxpoplplot=-1, intreval=False):
    
    if intreval:
        figrsize = [10, 10]
    else:
        figrsize = (gdat.sizeimag, gdat.sizeimag)
    figr, axis = plt.subplots(figsize=figrsize)
    
    if intreval:
        axis.set_position([0.5, 0.1, 0.6, 0.8])
        path = ''
    else:
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

        nameplot = '%s%s%s%s' % (strgplot, strgener, strgevtt, strgpopl)
   
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
    axis.axvline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axvline(-innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(innr, ls='--', alpha=gdat.alphbndr, color='black')
    axis.axhline(-innr, ls='--', alpha=gdat.alphbndr, color='black')


def retr_scat(gdat, axis, maps, thisindxener, thisindxevtt):

    draw_frambndr(gdat, axis)
    
    scat = axis.scatter(maps[thisindxener, :, thisindxevtt, 0], maps[thisindxener, :, thisindxevtt, 1], alpha=gdat.alphmaps, facecolor='black', s=5)

    return scat


def retr_imag(gdat, axis, maps, strg, strgcbar, thisindxener=None, thisindxevtt=-1, tdim=False, imag=None):
    
    vmin = getattr(gdat, 'minmscal' + strgcbar)
    vmax = getattr(gdat, 'maxmscal' + strgcbar)
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
        
        shapflat = list(maps.shape)
        shapflat[0] = gdat.numbpixlfull
        mapstemp = zeros(shapflat)
        mapstemp[gdat.indxpixlrofi, ...] = maps
        maps = mapstemp.reshape(shap).swapaxes(0, 1)
        
        #maps = maps.reshape(shap).swapaxes(0, 1)
    
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


def make_catllabl(gdat, strg, axis):
    
    # transdimensional elements
    if strg == 'post' and gdat.condcatl or strg == 'this':
        if strg == 'post':
            colr = 'black'
            labl = 'Condensed Model %s' % gdat.strgelem
        else:
            colr = 'b'
            labl = 'Sample Model %s' % gdat.strgelem
        axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                                                                            label=labl, marker='+', lw=gdat.mrkrlinewdth, color=colr)
        
    for q in gdat.indxrefr:
        if gdat.refrlgal[q] != None and gdat.refrbgal[q] != None:
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, \
                                                                                      label=gdat.legdrefrhits[q], marker='x', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                      label=gdat.legdrefrmiss[q], marker='s', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
    
    # fixed-dimensional objects
    if gdat.elemtype == 'lens':
        if strg == 'this':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='Model Source', marker='<', lw=gdat.mrkrlinewdth, color='b')
    
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='Model Host', marker='s', lw=gdat.mrkrlinewdth, color='b')
        if gdat.datatype == 'mock':
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='%s Source' % gdat.legdrefr[q], marker='>', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
        
            axis.scatter(gdat.anglfact * gdat.maxmgangdata * 5., gdat.anglfact * gdat.maxmgangdata * 5, s=50, alpha=gdat.alphelem, facecolor='none', \
                                                                                 label='%s Host' % gdat.legdrefr[q], marker='D', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
    
    temphand, temp = axis.get_legend_handles_labels()
    numblabl = len(temp)
    
    if numblabl == 4:
        numbcols = 2
    else:
        numbcols = 3
    axis.legend(bbox_to_anchor=[0.5, 1.15], loc='center', ncol=numbcols)
        

def supr_fram(gdat, gdatmodi, strg, axis, indxpoplplot=-1):
    
    strgmodl = retr_strgmodl(strg)
    
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    
    # associations with the reference elements
    if gdat.refrlgal != None and gdat.refrbgal != None: 
        for q in gdat.indxrefr:
            if gdat.refrnumbelem[q] > 0:
                if indxpoplplot == -1:
                    listindxpoplplot = gdat.indxrefr
                else:
                    listindxpoplplot = [indxpoplplot]
                for l in listindxpoplplot:
                    reframpl = getattr(gdat, 'refr' + gdat.namecompampl)
                    if reframpl == None:
                        mrkrsize = full(gdat.refrnumbelem[q], 5.)
                    else:
                        mrkrsize = retr_mrkrsize(gdat, reframpl[q][0, :])
                    lgal = copy(gdat.refrlgal[q][0, :])
                    bgal = copy(gdat.refrbgal[q][0, :])
                    numbelem = int(gdat.refrnumbelem[q])
                    if gdatmodi != None and numbtrap > 0:
                        ### hit
                        indx = gdatmodi.thisindxelemrefrasschits[q][l]
                        if indx.size > 0:
                            axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, label=gdat.legdrefrmiss, \
                                                                                                                      marker='x', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
                        ### missed
                        indx = gdatmodi.thisindxelemrefrasscmiss[q][l]
                    else:
                        indx = arange(lgal.size)
                    if indx.size > 0: 
                        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.alphelem, facecolor='none', \
                                                                                            label=gdat.legdrefrhits, marker='s', lw=gdat.mrkrlinewdth, color=gdat.listcolrrefr[q])
            
    ## host galaxy position
    if hostemistype != 'none':
        axis.scatter(gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalhost], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalhost], facecolor='none', \
                                                                    #alpha=gdat.alphelem, \
                                                                    alpha=0.7, \
                                                                    label=gdat.legdrefrhits, s=300, marker='D', lw=gdat.mrkrlinewdth, color='g')
    if lensmodltype != 'none':
        ## host galaxy Einstein radius
        axis.add_patch(plt.Circle((gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalhost], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalhost]), \
                                                                        gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbeinhost], \
                                                                        edgecolor='g', facecolor='none', lw=gdat.mrkrlinewdth))
        
        ## source galaxy position
        axis.scatter(gdat.anglfact * gdat.truefixp[gdat.trueindxfixplgalsour], gdat.anglfact * gdat.truefixp[gdat.trueindxfixpbgalsour], facecolor='none', \
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
                mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp[gdat.namecompampl][l]])
                lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['lgal'][l]]
                bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampcomp['bgal'][l]]
                axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.alphelem, label='Sample', marker='+', lw=gdat.mrkrlinewdth, color='b')

            ## source
            if lensmodltype != 'none':
                lgalsour = gdatmodi.thissampvarb[gdat.fittindxfixplgalsour]
                bgalsour = gdatmodi.thissampvarb[gdat.fittindxfixpbgalsour]
                axis.scatter(gdat.anglfact * lgalsour, gdat.anglfact * bgalsour, facecolor='none', \
                                                                      alpha=gdat.alphelem, \
                                                                      label='Model Source', s=300, marker='<', lw=gdat.mrkrlinewdth, color='b')
    
            if hostemistype != 'none':
                ## host
                lgalhost = gdatmodi.thissampvarb[gdat.fittindxfixplgalhost]
                bgalhost = gdatmodi.thissampvarb[gdat.fittindxfixpbgalhost]
                axis.scatter(gdat.anglfact * lgalhost, gdat.anglfact * bgalhost, facecolor='none', \
                                                                      alpha=gdat.alphelem, \
                                                                      label='Model Host', s=300, marker='s', lw=gdat.mrkrlinewdth, color='b')
                if lensmodltype != 'none':
                    beinhost = gdatmodi.thissampvarb[gdat.fittindxfixpbeinhost]
                    axis.add_patch(plt.Circle((gdat.anglfact * lgalhost, gdat.anglfact * bgalhost), gdat.anglfact * beinhost, edgecolor='b', facecolor='none', \
                                                                                                                                         lw=gdat.mrkrlinewdth, ls='--'))
                
    # temp
    if strg == 'post' and gdat.condcatl and gdat.fittnumbtrap > 0:
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
            
            # temp
            if attrtemp == 'meanelempop0':
                attrtemp = 'meanelem'
            if attrtemp == 'fluxdistsloppop0':
                attrtemp = 'fluxdistslop'

            if attrtemp in gdat.fittliststrgfeatodimtotl:
                # temp
                if attrtemp == 'deltllik':
                    continue

                for l in gdat.fittindxpopl:
                    attrprim = attrtemp + 'pop%d' % l
                    listtemp = []
                    arrytemp = zeros((gdat.numbsamptotl, gdat.fittmaxmnumbelem[l])) + nan
                    for n in range(len(valu)):
                        arrytemp[n, :len(valu[n][l])] = valu[n][l]
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


def initchro(gdat, gdatmodi, strg, name):

    if strg == 'next':
        gdatmodi.thischro[gdat.indxchro[name]] = gdat.functime()
    if strg == 'gene':
        gdatmodi.thischro[name] = gdat.functime()
    

def stopchro(gdat, gdatmodi, strg, name):
    
    if strg == 'next':
        gdatmodi.thischro[gdat.indxchro[name]] = gdat.functime() - gdatmodi.thischro[gdat.indxchro[name]]
    if strg == 'gene':
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
def retr_deflelem_jitt(deflelem, indxpixl, lgalgrid, bgalgrid, numbelemconc, lgalconc, bgalconc, defsconc, ascaconc, acutconc):
    
    for k in range(numbelemconc):
        deflelem[:] += retr_defl_jitt(indxpixl, lgalgrid, bgalgrid, lgalconc[k], bgalconc[k], defsconc[k], 0., 0., asca=ascaconc[k], acut=acutconc[k])


def retr_strgmodl(strg):

    if strg == 'true':
        strgmodl = strg
    else:
        strgmodl = 'fitt'
    
    return strgmodl


def proc_samp(gdat, gdatmodi, strg, raww=False, fast=False):

    initchro(gdat, gdatmodi, strg, 'proc')

    if gdat.verbtype > 1:
        print 'proc_samp()'
        print 'strg'
        print strg
        print

    if gdatmodi != None:
        gdatobjt = gdatmodi
    else:
        gdatobjt = gdat

    strgmodl = retr_strgmodl(strg)

    # common dictionary
    dicttemp = {}
           
    # temp
    numbtrap = getattr(gdat, strgmodl + 'numbtrap')
    if numbtrap > 0:
        indxfixpmeanelem = getattr(gdat, strgmodl + 'indxfixpmeanelem')
        indxpopl = getattr(gdat, strgmodl + 'indxpopl')
        if gdat.elemtype == 'lght':
            minmflux = getattr(gdat, strgmodl + 'minmflux')
        if gdat.elemtype == 'lens':
            minmdefs = getattr(gdat, strgmodl + 'minmdefs')
        if gdat.elemtype == 'clus':
            minmnobj = getattr(gdat, strgmodl + 'minmnobj')
        numbpopl = getattr(gdat, strgmodl + 'numbpopl')
        numbcomp = getattr(gdat, strgmodl + 'numbcomp')
        spectype = getattr(gdat, strgmodl + 'spectype')
        liststrgfeatdefa = getattr(gdat, strgmodl + 'liststrgfeatdefa')
        liststrgfeattotl = getattr(gdat, strgmodl + 'liststrgfeattotl')
        liststrgfeatodimtotl = getattr(gdat, strgmodl + 'liststrgfeatodimtotl')
        liststrgcomp = getattr(gdat, strgmodl + 'liststrgcomp')
        listscalcomp = getattr(gdat, strgmodl + 'listscalcomp')
        liststrgcomptotl = getattr(gdat, strgmodl + 'liststrgcomptotl')
        oaxitype = getattr(gdat, strgmodl + 'oaxitype')
    
    backtype = getattr(gdat, strgmodl + 'backtype')
    psfntype = getattr(gdat, strgmodl + 'psfntype')
    convdiff = getattr(gdat, strgmodl + 'convdiff')
    lensmodltype = getattr(gdat, strgmodl + 'lensmodltype')
    listnamediff = getattr(gdat, strgmodl + 'listnamediff')
    hostemistype = getattr(gdat, strgmodl + 'hostemistype')
    psfnevaltype = getattr(gdat, strgmodl + 'psfnevaltype')
    
    # grab the sample vector
    sampvarb = getattr(gdatobjt, strg + 'sampvarb')
    
    if psfnevaltype != 'none':
        psfp = sampvarb[getattr(gdat, strgmodl + 'indxfixppsfp')]

    bacp = sampvarb[getattr(gdat, strgmodl + 'indxfixpbacp')]
    
    if numbtrap > 0:
        indxelemfull = list(getattr(gdatobjt, strg + 'indxelemfull'))
        # temp -- this may slow down execution
        indxsampcomp = retr_indxsampcomp(gdat, indxelemfull, strgmodl)
        setattr(gdatobjt, strg + 'indxsampcomp', indxsampcomp)
        indxfixpnumbelem = getattr(gdat, strgmodl + 'indxfixpnumbelem')
    
        numbelem = sampvarb[indxfixpnumbelem].astype(int)
        
        dictelem = [[] for l in indxpopl]
        for l in indxpopl:
            dictelem[l] = dict()
            for strgfeat in liststrgfeatdefa:
                dictelem[l][strgfeat] = []
            for strgcomp in liststrgcomp[l]:
                dictelem[l][strgcomp] = sampvarb[indxsampcomp[strgcomp][l]]
        
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
            for l in indxpopl:
                dictelem[l]['spec'] = retr_spec(gdat, dictelem[l]['flux'], dictelem[l]['sind'], dictelem[l]['curv'], dictelem[l]['expc'], spectype=spectype[l])
        
        for strgfeat in gdat.liststrgfeatconc:
            if strgfeat == 'spec':
                dicttemp[strgfeat + 'conc'] = concatenate([dictelem[l][strgfeat] for l in indxpopl], axis=1)
            else:
                dicttemp[strgfeat + 'conc'] = concatenate([dictelem[l][strgfeat] for l in indxpopl])

        numbelemconc = dicttemp['lgalconc'].size
    else:
        numbelemconc = 0

    if gdat.verbtype > 1:
        if numbtrap > 0:
            for l in indxpopl:
                print 'l'
                print l
                for strgcomp in liststrgcomp[l]:
                    print strgcomp
                    print dictelem[l][strgcomp]
            print

    ### loglikelihood
    initchro(gdat, gdatmodi, strg, 'llik')

    if gdat.calcllik:
        
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
                beinhost = 0.
            else:
                beinhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpbeinhost')]
        if hostemistype != 'none':
            ellphost = sampvarb[getattr(gdat, strgmodl + 'indxfixpellphost')]
            anglhost = sampvarb[getattr(gdat, strgmodl + 'indxfixpanglhost')]
            serihost = sampvarb[getattr(gdat, strgmodl + 'indxfixpserihost')]
        if lensmodltype != 'none':
            if raww:
                sherextr = 0.
            else:
                sherextr = sampvarb[getattr(gdat, strgmodl + 'indxfixpsherextr')]
            sangextr = sampvarb[getattr(gdat, strgmodl + 'indxfixpsangextr')]
           
            initchro(gdat, gdatmodi, strg, 'deflzero')
            defl = zeros((gdat.numbpixl, 2))
            stopchro(gdat, gdatmodi, strg, 'deflzero')
            
            ## host halo deflection
            initchro(gdat, gdatmodi, strg, 'deflhost')
            if strg == 'next' and not gdatmodi.prophost:
                # retrieve state variable
                deflhost = gdatmodi.thisdeflhost
            else:
                deflhost = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost)
                setattr(gdatobjt, strg + 'deflhost', deflhost)
            defl += deflhost
            stopchro(gdat, gdatmodi, strg, 'deflhost')

            ## external shear
            initchro(gdat, gdatmodi, strg, 'deflextr')
            deflextr = retr_deflextr(gdat, sherextr, sangextr)
            defl += deflextr
            stopchro(gdat, gdatmodi, strg, 'deflextr')
        
        ## construct the PSF to be convolved with the image
        if psfnevaltype == 'conv' or psfnevaltype == 'full':
            initchro(gdat, gdatmodi, strg, 'psfnconv')
            if strg == 'next' and not gdatmodi.proppsfnconv:
                # retrieve state variable
                psfnconv = gdatmodi.thispsfnconv
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
                setattr(gdatobjt, strg + 'psfp', psfp)
                setattr(gdatobjt, strg + 'psfnconv', psfnconv)
            stopchro(gdat, gdatmodi, strg, 'psfnconv')
            
        if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and numbtrap > 0:
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
    
        # prepare the evaluation arrays
        if numbtrap > 0:
            initchro(gdat, gdatmodi, strg, 'pixlkern')
            if gdat.elemtype == 'lens':

                if strg == 'next' and gdat.propwithsing and gdat.pertmodleval:
                    # retrieve state variable
                    if gdatmodi.propfixp:
                        numbelemeval = 0
                        deflelem = gdatmodi.thisdeflelem
                    else:
                        numbelemeval = gdatmodi.numbelemeval
                        deflelem = copy(gdatmodi.thisdeflelem)
                        for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                            dicttemp[namecomp + 'eval'] = gdatmodi.dicttemp[namecomp + 'eval']
                else:
                    numbelemeval = numbelemconc
                    deflelem = zeros((gdat.numbpixl, 2))
                    for namecomp in getattr(gdat, strgmodl + 'liststrgcomptotl'):
                        dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']
            
            if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
                if strg == 'next' and gdat.pertmodleval and gdat.propwithsing and not gdatmodi.proppsfp:
                    # perturbative point source
                    if gdatmodi.propfixp:
                        numbelemeval = 0
                    if gdatmodi.propcomp or gdatmodi.proptran:
                        numbelemeval = gdatmodi.numbelemeval
                        if gdat.elemtype == 'lght':
                            dicttemp['sindeval'] = []
                            dicttemp['curveval'] = []
                            dicttemp['expceval'] = []
                        
                        for namecomp in gdat.fittliststrgcomp[gdatmodi.indxpoplmodi]:
                            dicttemp[namecomp + 'eval'] = gdatmodi.dicttemp[namecomp + 'eval']
                        
                        if gdat.elemtype == 'lght':
                            dicttemp['speceval'] = retr_spec(gdat, dicttemp['fluxeval'], dicttemp['sindeval'], dicttemp['curveval'], \
                                                                                                                    dicttemp['expceval'], spectype[gdatmodi.indxpoplmodi])
                else:
                    # all point sources
                    if gdat.verbtype > 1:
                        print 'Preparing all point sources...'
                    numbelemeval = numbelemconc
                    for namecomp in gdat.listnamefeateval:
                        dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']

            stopchro(gdat, gdatmodi, strg, 'pixlkern')
                
            if gdat.verbtype > 1:
                if numbelemeval > 0:
                    for namecomp in gdat.listnamefeateval:
                        print 'dicttemp[%seval]' % namecomp
                        print dicttemp[namecomp + 'eval']
                    print 'gdat.evalcirc'
                    print gdat.evalcirc
                    print
            
        # determine the pixels over which element kernels will be evaluated
        if gdat.evalcirc == 'full':
            if numbtrap > 0:
                listindxpixleval = gdat.listindxpixl[numbelemeval]
                listindxpixlevalconc = gdat.indxpixl
            indxpixleval = gdat.indxpixl
        if gdat.evalcirc == 'locl':
            if strg == 'next':
                if gdatmodi.evalllikpert:
                    listindxpixleval, listindxpixlevalconc = retr_indxpixlevalconc(gdat, dicttemp, gdat.evalcirc)
                    indxpixleval = listindxpixlevalconc
                else:
                    indxpixleval = gdat.indxpixl
            else:
                if numbtrap > 0:
                    listindxpixleval, listindxpixlevalconc = retr_indxpixlevalconc(gdat, dicttemp, gdat.evalcirc)
                indxpixleval = gdat.indxpixl
        if strg == 'next':
            gdatmodi.indxpixleval = indxpixleval
    
        if strg == 'next':
            # temp
            gdatmodi.indxenermodi = gdat.indxener
            gdatmodi.indxevttmodi = gdat.indxevtt
            gdatmodi.indxcubemodi = meshgrid(gdatmodi.indxenermodi, gdatmodi.indxpixleval, gdatmodi.indxevttmodi, indexing='ij')
            indxcube = gdatmodi.indxcubemodi
        else:
            indxcube = gdat.indxcube
        
        # element kernel evaluation
        if numbtrap > 0:
            initchro(gdat, gdatmodi, strg, 'kernelem')
            if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
                if strg == 'next' and gdat.propwithsing and gdat.pertmodleval and not gdatmodi.proppsfp:
                    if gdatmodi.propfixp:
                        sbrtpnts = gdatmodi.thissbrtpnts
                    if gdatmodi.propcomp or gdatmodi.proptran:
                        gdat.sbrtpntstemp[gdatmodi.indxcubemodi] = copy(gdatmodi.thissbrtpnts[gdatmodi.indxcubemodi])
                        sbrtpnts = gdat.sbrtpntstemp
                else:
                    sbrtpnts = zeros_like(gdat.expo)

                for k in range(numbelemeval):
                    if gdat.elemtype == 'lght':
                        varbevalextd = dicttemp['speceval'][:, k]
                    if gdat.elemtype == 'clus':
                        varbevalextd = dicttemp['nobjeval'][None, k]
                    if gdat.verbtype > 1:
                        print 'varbevalextd'
                        print varbevalextd
                    sbrtpnts[:, listindxpixleval[k], :] += retr_sbrtpnts(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], varbevalextd, \
                                                                                                                        psfnintp, oaxitype, listindxpixleval[k])
               
                # when the only background template is the data-PS residual, correct the PS template for numerical noise
                if backtype[0] == 'data':
                    print 'Correcting the PS sbrt...'
                    print 'sbrtpnts'
                    summgene(sbrtpnts)
                    sbrtpnts[where(sbrtpnts <= 1e-100)] = 1e-100
                    print 'sbrtpnts'
                    summgene(sbrtpnts)

                setattr(gdatobjt, strg + 'sbrtpnts', sbrtpnts)

                if gdat.diagmode:
                    if gdat.elemtype == 'lght' or gdat.elemtype == 'clus':
                        if amin(sbrtpnts) / mean(sbrtpnts) < -1e-10:
                            raise Exception('')
            
            if gdat.elemtype == 'lens':
                initchro(gdat, gdatmodi, strg, 'kernelem')
                
                if strg == 'next' and gdat.propwithsing and gdat.pertmodleval:
                    if gdatmodi.propfixp:
                        deflelem = gdatmodi.thisdeflelem
                    if gdatmodi.propcomp or gdatmodi.proptran:
                        deflelem = copy(gdatmodi.thisdeflelem)
                else:
                    deflelem = zeros((gdat.numbpixl, 2))

                for k in range(numbelemeval):
                    if gdat.variasca:
                        asca = dicttemp['ascaeval'][k]
                    else:
                        asca = gdat.ascaglob
                    if gdat.variacut:
                        acut = dicttemp['acuteval'][k]
                    else:
                        acut = gdat.acutglob
                    deflelem[listindxpixleval[k], :] += retr_defl_jitt(gdat.indxpixl, gdat.lgalgrid, gdat.bgalgrid, \
                                                 dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], dicttemp['defseval'][k], 0., 0., \
                                                                                                                asca=asca, acut=acut, indxpixleval=listindxpixleval[k])
                defl += deflelem
                setattr(gdatobjt, strg + 'deflelem', deflelem)
                stopchro(gdat, gdatmodi, strg, 'kernelem')
            stopchro(gdat, gdatmodi, strg, 'kernelem')

        sbrtdiff = dict()
        if lensmodltype != 'none':
            
            # lensed surface brightness
            initchro(gdat, gdatmodi, strg, 'sbrtlens')
            if gdat.numbener > 1:
                specsour = retr_spec(gdat, array([fluxsour]), array([sindsour]))
            else:
                specsour = array([[fluxsour]])
            sbrtdiff['lens'] = retr_sbrtsers(gdat, gdat.lgalgrid - defl[:, 0], gdat.bgalgrid - defl[:, 1], lgalsour, bgalsour, specsour, sizesour, ellpsour, anglsour)

            stopchro(gdat, gdatmodi, strg, 'sbrtlens')
            
            defl = defl.reshape((gdat.numbsidecart, gdat.numbsidecart, 2))
            setattr(gdatobjt, strg + 'defl', defl)
            
        ### background surface brightness
        numbback = getattr(gdat, strgmodl + 'numbback')
        indxback = getattr(gdat, strgmodl + 'indxback')
        sbrtbacknorm = getattr(gdat, strgmodl + 'sbrtbacknorm')
        indxbacpback = getattr(gdat, strgmodl + 'indxbacpback')
        indxtessback = getattr(gdat, strgmodl + 'indxtessback')
        # temp
        #smthback = getattr(gdat, strgmodl + 'smthback')
        specback = getattr(gdat, strgmodl + 'specback')
        unifback = getattr(gdat, strgmodl + 'unifback')
        sbrtback = empty((numbback, gdat.numbener, indxpixleval.size, gdat.numbevtt))
        
        # evaluate host galaxy surface brightness
        if hostemistype != 'none':
            initchro(gdat, gdatmodi, strg, 'sbrthost')
            if not strg == 'next' or gdatmodi.prophost:
                if gdat.verbtype > 1:
                    print 'Host galaxy surface brightness evaluation...'
                if gdat.numbener > 1:
                    spechost = retr_spec(gdat, array([fluxhost]), array([sindhost]))
                else:
                    spechost = array([[fluxhost]])
                sbrtdiff['host'] = retr_sbrtsers(gdat, gdat.lgalgrid, gdat.bgalgrid, lgalhost, bgalhost, spechost, sizehost, ellphost, anglhost, serihost)
                setattr(gdatobjt, strg + 'sbrthost', sbrtdiff['host'])
            if gdat.verbtype > 1:
                print 'sbrtdiff[host]'
                summgene(sbrtdiff['host'])
            stopchro(gdat, gdatmodi, strg, 'sbrthost')
        
        # construct the model diffuse surface brightness
        # convolve the model surface brightness with the PSF
        if gdat.verbtype > 1:
            print 'PSF convolution of diffuse components...'

        initchro(gdat, gdatmodi, strg, 'sbrtdiffconv')
        sbrtdiffconv = dict()
        for k, name in enumerate(listnamediff):
            
            if gdat.verbtype > 1:
                print name
           
            booltemp = False
            if strg == 'next':
                if name.startswith('back'):
                    if not gdatmodi.proppsfp or unifback[int(name[4:])]:
                        booltemp = True
                elif name == 'host':
                    if not gdatmodi.prophost and not gdatmodi.proppsfp:
                        booltemp = True
            
            if booltemp:
                sbrtdiffconv[name] = getattr(gdatmodi, 'thissbrt' + name + 'conv')
                if gdat.verbtype > 1:
                    print 'Retrieving diffuse component from the state vector...'
            else:
                
                if name.startswith('back'):
                    indxbacktemp = int(name[4:])
                    if gdat.pixltype == 'heal' and (psfnevaltype == 'full' or psfnevaltype == 'conv') and not unifback[indxbacktemp]:
                        sbrtdiff[name] = getattr(gdat, strgmodl + 'sbrtbackhealfull')[indxbacktemp, :, :, :]
                    else:
                        sbrtdiff[name] = sbrtbacknorm[indxbacktemp, :, :, :]
                
                if convdiff[k] and (psfnevaltype == 'full' or psfnevaltype == 'conv') and gdat.lensmodltype != 'none':
                    #sbrtdiffconv[name] = empty((gdat.numbener, indxpixleval.size, gdat.numbevtt))
                    sbrtdiffconv[name] = empty_like(gdat.expo)
                    # temp -- change this to modified energy and event class bins
                    if gdat.pixltype == 'heal':
                        psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
                        fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            if gdat.pixltype == 'cart':
                                sbrtdiffconv[name][i, :, m] = convolve_fft(sbrtdiff[name][i, :, m].reshape((gdat.numbsidecart, gdat.numbsidecart)), psfnconv[m][i]).flatten()
                                indx = where(sbrtdiffconv[name][i, :, m] < 1e-6)
                                sbrtdiffconv[name][i, indx, m] = 1e-6
                            if gdat.pixltype == 'heal':
                                sbrtdiffconv[name][i, :, m] = hp.smoothing(sbrtdiff[name][i, :, m], fwhm=fwhm[i, m])[gdat.indxpixlrofi]
                                sbrtdiffconv[name][i, :, m][where(sbrtdiffconv[name][i, :, m] <= 1e-50)] = 1e-50
                    if gdat.verbtype > 1:
                        print 'Convolving...'
                else:
                    sbrtdiffconv[name] = sbrtdiff[name]
                    if gdat.verbtype > 1:
                        print 'Falling back on the unconvolved components...'
                setattr(gdatobjt, strg + 'sbrt' + name + 'conv', sbrtdiffconv[name])
                if gdat.verbtype > 1:
                    print 'sbrtdiff[name]'
                    summgene(sbrtdiff[name])
            if gdat.verbtype > 1:
                print 'sbrtdiffconv[name]'
                summgene(sbrtdiffconv[name])
    
        stopchro(gdat, gdatmodi, strg, 'sbrtdiffconv')
        
        initchro(gdat, gdatmodi, strg, 'sbrtmodl')
        ## initialize with the must-have background
        if gdat.verbtype > 1:
            print 'Summing up the model emission...'
        if specback[0]:
            sbrtmodl = sbrtdiffconv['back0000'][indxcube] * bacp[indxbacpback[0]]
        else:
            sbrtmodl = sbrtdiffconv['back0000'][indxcube] * bacp[indxbacpback[0]][:, None, None]
        
        ## add the other diffuse emission components
        for name in listnamediff:
            if name == 'back0000':
                continue
            elif name.startswith('back'):
                indxbacktemp = int(name[4:])
                if specback[indxbacktemp]:
                    sbrtmodl += sbrtdiffconv[name][indxcube] * bacp[indxbacpback[indxbacktemp]]
                else:
                    sbrtmodl += sbrtdiffconv[name][indxcube] * bacp[indxbacpback[indxbacktemp]][:, None, None]
            else:
                sbrtmodl += sbrtdiffconv[name][indxcube]
            if gdat.verbtype > 1:
                print name
                print 'sbrtdiffconv[name][indxcube]'
                summgene(sbrtdiffconv[name][indxcube])
                print 'sbrtmodl'
                summgene(sbrtmodl)
                
        ## add point source emission
        if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and numbtrap > 0:
            sbrtmodl += sbrtpnts[indxcube]
        stopchro(gdat, gdatmodi, strg, 'sbrtmodl')
            
        ### count map
        initchro(gdat, gdatmodi, strg, 'expo')
        cntpmodl = retr_cntp(gdat, sbrtmodl, indxpixlmean=indxpixleval)
        
        setattr(gdatobjt, strg + 'cntpmodl', cntpmodl)
        stopchro(gdat, gdatmodi, strg, 'expo')

        # mock data specific
        if strg == 'true':
            
            if raww:
                strgvarb = 'truecntpmodlraww'
            else:
                strgvarb = 'cntpdata'
            
            # generate count data
            cntstemp = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
            for i in gdat.indxener:
                for j in gdat.indxpixl:
                    for m in gdat.indxevtt:
                        cntstemp[i, j, m] = poisson(cntpmodl[i, j, m])
            setattr(gdat, strgvarb, cntstemp)
            
            if raww:
                return
            else:
                retr_datatick(gdat)
       
            proc_cntpdata(gdat)
        
        # diagnostics
        if gdat.diagmode:
            if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and numbtrap > 0:
                if amin(sbrtpnts) < -mean(sbrtpnts) * 1e-10:
                    print 'sbrtpnts'
                    summgene(sbrtpnts)
                    raise Exception('Point source spectral surface brightness is not positive-definite.')

            for varb in [sbrtmodl, cntpmodl]:
                frac = amin(varb) / mean(varb)
                if frac < -1e-10:
                    print 'frac'
                    print frac
                    summgene(varb)
                    raise Exception('Total model spectral surface brightness is not positive-definite.')

        ### log-likelihood
        initchro(gdat, gdatmodi, strg, 'llikcalc')
        #lliktemp = retr_llik_mult(gdat, cntpmodl)
        
        lliktemp = retr_llik_depr(gdat, cntpmodl, indxpixleval)
        
        if gdat.verbtype > 1:
            print 'cntpmodl'
            summgene(cntpmodl)
            print 'indxpixleval'
            summgene(indxpixleval)
            print 'lliktemp'
            summgene(lliktemp)
            print

        stopchro(gdat, gdatmodi, strg, 'llikcalc')
        
        if gdat.diagmode:
            if not isfinite(lliktemp).all():
                raise Exception('Likelihood is not finite.')
    
        if strg == 'next':
            gdatmodi.thisdeltlliktotl = sum(lliktemp - gdatmodi.thisllik[:, indxpixleval, :])
            lliktotl = gdatmodi.thisdeltlliktotl + gdatmodi.thislliktotl
        else:
            lliktotl = sum(lliktemp)
        
        setattr(gdatobjt, strg + 'llik', lliktemp)
    else:
        lliktotl = 0.
    
    setattr(gdatobjt, strg + 'lliktotl', lliktotl) 
    
    stopchro(gdat, gdatmodi, strg, 'llik')

    # log-prior
    initchro(gdat, gdatmodi, strg, 'lpri')
    lpri = zeros(gdat.numblpri)
    if numbtrap > 0:
        
        liststrgfeatprio = getattr(gdat, strgmodl + 'liststrgfeatprio')
        liststrgpdfnprio = getattr(gdat, strgmodl + 'liststrgpdfnprio')
    
        meanelem = sampvarb[indxfixpmeanelem]
        
        if gdat.penalpridiff:
            
            if strg == 'next' and gdatmodi.prophypr:
                lpri[1] = gdatmodi.thislpridiff
            else:
                sbrtdatapnts = gdat.sbrtdata - sbrtpnts
                if gdat.pixltype == 'heal':
                    raise Exception('')
                if gdat.pixltype == 'cart':
                    psecodimdatapnts = empty((gdat.numbener, gdat.numbsidecart / 2, gdat.numbevtt))
                    psfn = retr_psfn(gdat, psfp, gdat.indxener, gdat.binsangl, psfntype, gdat.binsoaxi, oaxitype)
                    fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
                    sigm = fwhm / 2.355
                    psecodimdatapntsprio = exp(-2. * gdat.meanmpolodim[None, :, None] / (1. / sigm[:, None, :]))
                    lpridiff = 0.
                    for i in gdat.indxener:
                        for m in gdat.indxevtt:
                            psecdatapnts = retr_psec(gdat, sbrtdatapnts[i, :, m].reshape((gdat.numbsidecart, gdat.numbsidecart)))
                            psecodimdatapnts[i, :, m] = retr_psecodim(gdat, psecdatapnts)
                            psecodimdatapnts[i, :, m] /= psecodimdatapnts[i, 0, m]
                            lpridiff += -0.5 * sum((psecodimdatapnts[i, :, m] - psecodimdatapntsprio[i, :, m])**2)
                            setattr(gdatobjt, strg + 'psecodimdatapntsene%devt%d' % (i, m), psecodimdatapnts[i, :, m])
                            setattr(gdatobjt, strg + 'psecodimdatapntsprioene%devt%d'% (i, m), psecodimdatapntsprio[i, :, m])
                lpri[1] = lpridiff 
                if strg == 'this' or strg == 'next':
                    setattr(gdatobjt, strg + 'lpridiff', lpridiff)
        
        for l in gdat.fittindxpopl:

            lpri[0] -= 0.5 * gdat.priofactdoff * numbcomp[l] * numbelem[l]
            lpri[2+0*numbpopl+l] = retr_probpois(numbelem[l], meanelem[l])
            
            for k, (strgfeat, pdfnfeat) in enumerate(zip(liststrgfeatprio[l], liststrgpdfnprio[l])):
                
                if False and pdfnfeat == 'tmpl':

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
                
                if False and pdfnfeat == 'self':
                    minmfeat = getattr(gdat, 'minm' + strgfeat)
                    maxmfeat = getattr(gdat, 'maxm' + strgfeat)
                    # temp -- this may be sped up a bit
                    lpri[2+(k+1)*numbpopl+l] = numbelem[l] * log(1. / (maxmfeat - minmfeat))
                
                if False:
                    if strgpdfn == 'disc':
                        indxfixpbgaldistscal = getattr(gdat, strgmodl + 'indxfixpbgaldistscalpop%d' % l)
                        lpri[2+2*numbpopl+l] = sum(log(pdfn_dexp(dictelem[l]['bgal'], gdat.maxmgang, sampvarb[indxfixpbgaldistscal]))) 
                    elif strgpdfn == 'exposcal':
                        gang = retr_gang(dictelem[l]['lgal'], dictelem[l]['bgal'])
                        indxfixpgangdistscal = getattr(gdat, strgmodl + 'indxfixpgangdistscalpop%d' % l)
                        lpri[2+numbpopl+l] = sum(log(pdfn_expo(gang, gdat.maxmgang, sampvarb[indxfixpgangdistscal]))) 
                        lpri[2+2*numbpopl+l] = -numbelem[l] * log(2. * pi) 
                    elif strgpdfn == 'tmpl':
                        lpri[2+numbpopl+l] = sum(lpdfspatprioobjt(dictelem[l]['bgal'], dictelem[l]['lgal'], grid=False))
               
                # temp -- this should not be here
                if pdfnfeat == 'powrslop':
                    lpri[2+(k+1)*numbpopl+l] = retr_lpripowrdist(gdat, gdatmodi, strgmodl, dictelem[l][strgfeat], strgfeat, sampvarb, l)
                #if pdfnfeat == 'gausmeanstdv':
                #    lpri[2+(k+1)*numbpopl+l] = retr_lprigausdist(gdat, gdatmodi, strgmodl, dictelem[l][strgfeat], strgfeat, sampvarb, l)
            
        if strg == 'this':
            gdatmodi.thislpripena = lpri[0]
        
        # temp
        if strg == 'next':
            gdatmodi.thisdeltlpri = 0.
            gdatmodi.thislpau = zeros(gdat.numblpau)
            # temp -- generalize this to explicit prior evaluations in the parameters 
            if gdatmodi.prophypr:
                gdatmodi.thisdeltlpri = sum(lpri) - gdatmodi.thislpritotl
            else:
                gdatmodi.thisdeltlpri = 0.
        # temp
        if False and strg == 'next' and gdatmodi.proptran:
            
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
            
            if gdatmodi.proptran:
                gdatmodi.thisdeltlpri = gdatmodi.thislpautotl + sum(lpri)
            
    
    lpritotl = sum(lpri)
    
    setattr(gdatobjt, strg + 'lpritotl', lpritotl) 
    setattr(gdatobjt, strg + 'lpri', lpri)
    
    stopchro(gdat, gdatmodi, strg, 'lpri')
    
    if not gdat.calcllik:
        return
    else:
        lpostotl = lpritotl + lliktotl
        setattr(gdatobjt, strg + 'lpostotl', lpostotl) 
    
    if strg == 'next':
        setattr(gdatmodi, 'thislpriprop', lpri)
    
    stopchro(gdat, gdatmodi, strg, 'proc')
    
    if fast:
        return
    
    ## tertiary variables that are not needed for evolving the chain
    if strg != 'next':
       
        ## load necessary variables
        numbdeflsubhplot = getattr(gdat, strgmodl + 'numbdeflsubhplot')
        numbdeflsingplot = getattr(gdat, strgmodl + 'numbdeflsingplot')
        liststrgfeatcorr = getattr(gdat, strgmodl + 'liststrgfeatcorr')
        liststrgfeatcorrtotl = getattr(gdat, strgmodl + 'liststrgfeatcorrtotl')
        numbback = getattr(gdat, strgmodl + 'numbback')
      
        ## set previously calculated tertiary variables
        setattr(gdatobjt, strg + 'bacp', bacp)
        if psfnevaltype != 'none':
            setattr(gdatobjt, strg + 'psfp', psfp)
        
        ## derived variables
        ## residual count map 
        cntpresi = gdat.cntpdata - cntpmodl
        setattr(gdatobjt, strg + 'cntpresi', cntpresi)
        
        ### spatial averages
        #### background components
        sbrtbackmean = []
        for c in indxback:
            sbrtbackmean.append(retr_spatmean(gdat, sbrtback[c]))
            setattr(gdatobjt, strg + 'sbrtback%04dmean' % c, sbrtbackmean[c])
        #### total model
        numblablsbrt = getattr(gdat, strgmodl + 'numblablsbrt')
        if numblablsbrt > 1:
            sbrtmodlmean = retr_spatmean(gdat, sbrtmodl)
            setattr(gdatobjt, strg + 'sbrtmodlmean', sbrtmodlmean)
    
        ### count maps
        if (gdat.elemtype == 'lght' or gdat.elemtype == 'clus') and numbtrap > 0:
            if oaxitype:
                setattr(gdatobjt, strg + 'factoaxi', factoaxi)
           
            setattr(gdatobjt, strg + 'psfn', psfn)

            ### PSF FWHM
            fwhm = 2. * retr_psfnwdth(gdat, psfn, 0.5)
            setattr(gdatobjt, strg + 'fwhm', fwhm)
            
	    	### mean PS surface brightness
            sbrtpntsmean = retr_spatmean(gdat, sbrtpnts)
            setattr(gdatobjt, strg + 'sbrtpntsmean', sbrtpntsmean)

            if gdat.calcerrr and gdat.fittnumbtrap > 0:
                sbrtpntsfull = zeros_like(gdat.expo)
                numbelemeval = numbelemconc
                for namecomp in gdat.listnamefeateval:
                    dicttemp[namecomp + 'eval'] = dicttemp[namecomp + 'conc']
                for k in range(numbelemeval):
                    if gdat.elemtype == 'lght':
                        varbevalextd = dicttemp['speceval'][:, k]
                    if gdat.elemtype == 'clus':
                        varbevalextd = dicttemp['nobjeval'][None, k]
                    sbrtpntsfull[:, :, :] += retr_sbrtpnts(gdat, dicttemp['lgaleval'][k], dicttemp['bgaleval'][k], varbevalextd, psfnintp, oaxitype, gdat.indxpixl)
                cntppnts = retr_cntp(gdat, sbrtpntsfull)
                cntppntsfull = retr_cntp(gdat, sbrtpntsfull)
                cntperrr = cntppnts - cntppntsfull
                setattr(gdatobjt, strg + 'cntperrr', cntperrr)
                if amax(fabs(cntperrr / cntppntsfull)) > 0.1:
                    raise Exception('Approximation error in calculating the PS surface brightness is above the tolerance level.')

            #fluxbrgt, fluxbrgtassc = retr_fluxbrgt(gdat, dicttemp['lgalconc'], dicttemp['bgalconc'], concatenate(dicttemp['flux']))
            #setattr(gdatobjt, strg + 'fluxbrgt', fluxbrgt)
            #setattr(gdatobjt, strg + 'fluxbrgtassc', fluxbrgtassc)
    
        #### background 
        if gdat.correxpo:
            cntpback = empty((numbback, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
            cntpbacktotl = zeros_like(gdat.expo)
            for c in indxback:
                for name in listnamediff:
                    if name == 'back%04d' % c:
                        if specback[c]:
                            sbrt = sbrtdiffconv[name] * bacp[indxbacpback[c]]
                        else:
                            sbrt = sbrtdiffconv[name] * bacp[indxbacpback[c]][:, None, None]
                cntpback[c, ...] = retr_cntp(gdat, sbrt)
                cntpbacktotl += cntpback[c, ...] 
            setattr(gdatobjt, strg + 'cntpback', cntpback)
            setattr(gdatobjt, strg + 'cntpbacktotl', cntpbacktotl)
    
        cntpdiffconv = dict()
        cntpdiff = dict()
        for k, name in enumerate(listnamediff):
            
            indxpixltemp = gdat.indxpixl
            if name.startswith('back'):
                indxbacktemp = int(name[4:])
                if specback[indxbacktemp]:
                    fact = bacp[indxbacpback[indxbacktemp]]
                else:
                    fact = bacp[indxbacpback[indxbacktemp]][:, None, None]
                if gdat.pixltype == 'heal' and not unifback[indxbacktemp] and (psfnevaltype == 'full' or psfnevaltype == 'conv'):
                    indxpixltemp = gdat.indxpixlrofi
            else:
                fact = 1.
            cntpdiff[name] = retr_cntp(gdat, sbrtdiff[name][:, indxpixltemp, :] * fact)
            setattr(gdatobjt, strg + 'cntp' + name, cntpdiff[name])
            
            sbrtmean = retr_spatmean(gdat, sbrtdiff[name][:, indxpixltemp, :] * fact)
            setattr(gdatobjt, strg + 'sbrt' + name + 'mean', sbrtmean)
            
            if psfnevaltype == 'conv' and convdiff[k]:
                cntpdiffconv[name] = retr_cntp(gdat, sbrtdiffconv[name] * fact)
                setattr(gdatobjt, strg + 'cntp' + name + 'conv', cntpdiffconv[name])
                
        if lensmodltype != 'none':
            
            if strgmodl == 'true':
                s2nr = cntpdiffconv['lens'] / sqrt(cntpmodl)
                setattr(gdatobjt, strg + 's2nr', s2nr)
            
            masshostbein = array([gdat.massfrombein * beinhost**2])
            setattr(gdatobjt, strg + 'masshostbein', masshostbein)
            if numbtrap > 0:
                ### sort with respect to deflection at scale radius
                if numbelemconc > 0:
                    indxelemsortbrgt = argsort(dicttemp[gdat.namefeatsort + 'conc'])[::-1]
                    for strgcomp in liststrgcomptotl:
                        dicttemp[strgcomp + 'sort'] = dicttemp[strgcomp + 'conc'][indxelemsortbrgt][:numbelemconc]

            deflsing = zeros((gdat.numbpixl, 2, numbdeflsingplot))
            numbdeflsing = min(numbdeflsubhplot, numbelemconc) + 2
            if numbelemconc > 0:
                numbdeflsing += 1
            for k in range(numbdeflsing):
                if k == 0:
                    deflhost = retr_defl(gdat, lgalhost, bgalhost, beinhost, ellphost, anglhost)
                    deflsing[:, :, k] = deflhost
                elif k == 1:
                    deflsing[:, :, k] = deflextr
                elif k == 2:
                    deflsing[:, :, k] = defl.reshape((gdat.numbpixl, 2)) - deflextr - deflhost
                else:
                    
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
    
        if lensmodltype != 'none':
            sbrtlensmean = mean(sum(sbrtdiff['lens'], 2), 0)
            deflmgtd = sqrt(sum(defl**2, axis=2))
            cntplensgrad = empty((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt, 2))
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    cntplensgrad[i, :, :, m, :] = retr_gradmaps(gdat, cntpdiffconv['lens'][i, :, m]) * gdat.sizepixl
            cntplensgradmgtd = sqrt(sum(cntplensgrad**2, axis=4))
            cntplensgrad *= gdat.sizepixl
            indx = where(fabs(cntplensgrad) > 1. * gdat.sizepixl)
            cntplensgrad[indx] = sign(cntplensgrad[indx]) * 1. * gdat.sizepixl
            setattr(gdatobjt, strg + 'deflmgtd', deflmgtd)
            setattr(gdatobjt, strg + 'cntplensgrad', cntplensgrad)
            setattr(gdatobjt, strg + 'cntplensgradmgtd', cntplensgradmgtd)

        ## element related
        if numbtrap > 0:
            if gdat.datatype == 'mock' and strg == 'true':
                gdat.refrlgal = []
                gdat.refrbgal = []
                for l in gdat.trueindxpopl:
                    gdat.refrlgal.append(tile(dictelem[l]['lgal'], [3] + list(ones(dictelem[l]['lgal'].ndim, dtype=int))))
                    gdat.refrbgal.append(tile(dictelem[l]['bgal'], [3] + list(ones(dictelem[l]['bgal'].ndim, dtype=int))))
                gdat.refrnumbelem = gdat.truenumbelem
    
            # correlate the fitting model elements with the reference elements
            if gdat.refrinfo and gdat.refrlgal != None and gdat.refrbgal != None and not (strg == 'true' and gdat.datatype == 'mock'):
                indxelemrefrasschits = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                indxelemfittasschits = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        indxelemfittmatr = empty((gdat.refrnumbelem[q], numbelem[l]), dtype=int)
                        indxelemrefrmatr = empty((gdat.refrnumbelem[q], numbelem[l]), dtype=int)
                        matrdist = empty((gdat.refrnumbelem[q], numbelem[l]))
                        for k in range(numbelem[l]):
                            # construct a matrix of angular distances between reference and fitting elements
                            matrdist[:, k] = retr_angldist(gdat, gdat.refrlgal[q][0, :], gdat.refrbgal[q][0, :], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k])
                            indxelemrefrmatr[:, k] = arange(gdat.refrnumbelem[q])
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
                            if indxelemrefrmatr[c] in indxelemrefrasschits[q][l] or indxelemfittmatr[c] in indxelemfittasschits[q][l]:
                                continue
                            indxelemrefrasschits[q][l].append(indxelemrefrmatr[c])
                            indxelemfittasschits[q][l].append(indxelemfittmatr[c])
                        
                        indxelemrefrasschits[q][l] = array(indxelemrefrasschits[q][l])
                        indxelemfittasschits[q][l] = array(indxelemfittasschits[q][l])
                setattr(gdatobjt, strg + 'indxelemrefrasschits', indxelemrefrasschits)
                setattr(gdatobjt, strg + 'indxelemfittasschits', indxelemfittasschits)

                indxelemrefrasscmiss = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                indxelemfittasscfals = [[[] for l in gdat.fittindxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        # indices of the reference elements not associated with the fitting model elements
                        if gdat.refrnumbelem[q] > 0:
                            indxelemrefrasscmiss[q][l] = setdiff1d(arange(gdat.refrnumbelem[q]), indxelemrefrasschits[q][l])
                        # indices of the fitting model elements not associated with the reference elements
                        if numbelem[l] > 0:
                            indxelemfittasscfals[q][l] = setdiff1d(arange(numbelem[l]), indxelemfittasschits[q][l])
                setattr(gdatobjt, strg + 'indxelemrefrasscmiss', indxelemrefrasscmiss)
                setattr(gdatobjt, strg + 'indxelemfittasscfals', indxelemfittasscfals)
                
                for l in gdat.fittindxpopl:
                    for q in gdat.indxrefr:
                        # collect the associated reference element feature for each fitting element 
                        for namefeatrefr in gdat.listnamefeatrefronly[q][l]:
                            name = namefeatrefr + gdat.listnamerefr[q]
                            refrfeat = getattr(gdat, 'refr' + namefeatrefr)
                            dictelem[l][strgfeat + gdat.listnamerefr[q]] = zeros(numbelem[l])
                            if len(indxelemrefrasschits[q][l]) > 0 and refrfeat != None:
                                dictelem[l][strgfeat + gdat.listnamerefr[q]][indxelemfittasscfals[q][l]] = refrfeat[q][0, indxelemrefrasschits[q][l]]
                        
            ### derived quantities
            for l in indxpopl:
                #### radial and angular coordinates
                dictelem[l]['gang'] = retr_gang(dictelem[l]['lgal'], dictelem[l]['bgal'])
                dictelem[l]['aang'] = retr_aang(dictelem[l]['lgal'], dictelem[l]['bgal'])
               
                if gdat.elemtype == 'lght':
                    #### number of expected counts
                    dictelem[l]['cnts'] = retr_cntspnts(gdat, dictelem[l]['lgal'], dictelem[l]['bgal'], dictelem[l]['spec'])
                
            #### delta log-likelihood
            for l in indxpopl:
                dictelem[l]['deltllik'] = zeros(numbelem[l])
            if gdat.calcllik and not (strgmodl == 'true' and gdat.checprio): 
                if gdat.verbtype > 1:
                    print
                    print 'Calculating log-likelihood differences when removing elements from the model.'
                gdatmoditemp = tdpy.util.gdatstrt()
                for l in indxpopl:
                    if gdat.verbtype > 1:
                        print 'l'
                        print l
                    for k in range(numbelem[l]):
                        prep_gdatmodi(gdat, gdatmoditemp, gdatobjt, strg)
                        retr_thisindxprop(gdat, gdatmoditemp, deth=True, thisindxpopl=l)
                        retr_prop(gdat, gdatmoditemp, thisindxelem=k)
                        proc_samp(gdat, gdatmoditemp, 'next')
                        dictelem[l]['deltllik'][k] = lliktotl - gdatmoditemp.nextlliktotl
                        if gdat.verbtype > 1:
                            print 'k'
                            print k
                            print 'lliktotl'
                            print lliktotl
                            print 'gdatmoditemp.nextlliktotl'
                            print gdatmoditemp.nextlliktotl
                            print
                    if gdat.verbtype > 1:
                        print
                if gdat.verbtype > 1:
                    print 'deltllik calculation ended.'
                    print

            if gdat.elemtype == 'lght':
                #### spectra
                for l in indxpopl:
                    dictelem[l]['specplot'] = retr_spec(gdat, dictelem[l]['flux'], dictelem[l]['sind'], dictelem[l]['curv'], dictelem[l]['expc'], spectype[l], plot=True)
                    
            if gdat.elemtype == 'lens':
                
                #### distance to the source
                for l in range(numbpopl):
                    dictelem[l]['diss'] = retr_angldist(gdat, dictelem[l]['lgal'],  dictelem[l]['bgal'], lgalsour, bgalsour)
                
                for l in indxpopl:
                    dictelem[l]['deflprof'] = empty((gdat.numbangl, numbelem[l]))
                    dictelem[l]['mcut'] = empty(numbelem[l])
                    dictelem[l]['rele'] = empty(numbelem[l])
                    dictelem[l]['reln'] = empty(numbelem[l])
                    dictelem[l]['reld'] = empty(numbelem[l])
                    dictelem[l]['relc'] = empty(numbelem[l])
                    deflsubh = zeros((gdat.numbpixl, 2, numbelem[l]))
                    for k in arange(numbelem[l]):
                        
                        if gdat.variasca:
                            asca = dictelem[l]['asca'][k]
                        else:
                            asca = gdat.ascaglob

                        if gdat.variacut:
                            acut = dictelem[l]['acut'][k]
                        else:
                            acut = gdat.acutglob
                        
                        #### deflection profiles
                        dictelem[l]['deflprof'][:, k] = retr_deflcutf(gdat.meanangl, dictelem[l]['defs'][k], asca, acut)
             
                        ### truncated mass 
                        dictelem[l]['mcut'][k] = retr_mcut(gdat, dictelem[l]['defs'][k], asca, acut)
                        
                        #### dot product with the source flux gradient
                        # temp -- weigh the energy and PSF bins
                        dictelem[l]['rele'][k] = retr_rele(gdat, cntpdiffconv['lens'][0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], 1., asca, acut, gdat.indxpixl)
                        dictelem[l]['reln'][k] = dictelem[l]['rele'][k] / dictelem[l]['defs'][k] * gdat.sizepixl
                        dictelem[l]['reld'][k] = retr_rele(gdat, gdat.cntpdata[0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], 1., asca, acut, gdat.indxpixl)
                        dictelem[l]['relc'][k] = -retr_rele(gdat, cntpdiffconv['lens'][0, :, 0], dictelem[l]['lgal'][k], dictelem[l]['bgal'][k], 1., \
                                                                                                                            asca, acut, gdat.indxpixl, absv=False)
                   
            if strgmodl == 'true':
                for namesele in gdat.listnamesele:
                    indx = [[] for l in indxpopl]
                    for l in indxpopl:
                        if namesele == 'pars':
                            indx[l] = where(dictelem[l]['deltllik'] > 0.5 * numbcomp[l])[0]
                        if namesele == 'nrel':
                            indx[l] = where(dictelem[l]['reln'] > 0.3)[0]
                        for strgfeat in gdat.trueliststrgfeat[l]:
                            if strgfeat in gdat.listnamefeatsele:
                                dictelem[l][strgfeat + namesele] = zeros(numbelem[l])
                                dictelem[l][strgfeat + namesele][indx[l]] = dictelem[l][strgfeat][indx[l]]
                    setattr(gdat, 'trueindxelem' + namesele, indx)
                
            ### distribution of element parameters and features
            #### calculate the model filter
            listindxfittelemfilt = []
            for namemodlfilt in gdat.listnamemodlfilt:
                indxfittelemfilt = [[] for l in gdat.fittindxpopl]
                if namemodlfilt == '':
                    for l in gdat.fittindxpopl:
                        indxfittelemfilt[l] = where((fabs(dictelem[l]['lgal']) < gdat.maxmgangdata) & (fabs(dictelem[l]['bgal']) < gdat.maxmgangdata))[0]
                listindxfittelemfilt.append(indxfittelemfilt)
            
            #### one dimensional
            for l in indxpopl:
                for strgfeat in liststrgfeat[l]:
                    if strgfeat == 'spec':
                        temp = zeros((gdat.numbbinsplot, gdat.numbener))
                    else:
                        temp = zeros(gdat.numbbinsplot)
                    dictelem[l]['hist' + strgfeat] = temp
                    if strgfeat == 'specplot' or strgfeat == 'deflprof':
                        continue
                    elif strgfeat == 'spec':
                        for i in gdat.indxener:
                            dictelem[l]['hist' + strgfeat][:, i] = histogram(dictelem[l]['spec'][i, listindxfittelemfilt[0][l]], gdat.binsspec)[0]
                    elif strgfeat == 'cnts':
                        dictelem[l]['hist' + strgfeat] = histogram(dictelem[l]['cnts'][listindxfittelemfilt[0][l]], gdat.binscnts)[0]
                    elif not (strgfeat == 'curv' and spectype[l] != 'curv' or strgfeat == 'expc' and spectype[l] != 'expc'):
                        bins = getattr(gdat, 'bins' + strgfeat)
                        if strgfeat[:-4] in gdat.listnamefeatsele and strgmodl == 'true':
                            indx = intersect1d(getattr(gdat, strgmodl + 'indxelem' + strgfeat[-4:])[l], listindxfittelemfilt[0][l])
                            if indx.size > 0:
                                dictelem[l]['hist' + strgfeat] = histogram(dictelem[l][strgfeat[:-4]][indx], bins)[0]
                        else:
                            if len(dictelem[l][strgfeat]) > 0:
                                # temp
                                try:
                                    dictelem[l]['hist' + strgfeat] = histogram(dictelem[l][strgfeat][listindxfittelemfilt[0][l]], bins)[0]
                                except:
                                    print 'hey'
                                    print 'histograming failed'
                                    print 'strgfeat'
                                    print strgfeat
                                    print 'dictelem[l][strgfeat]'
                                    print dictelem[l][strgfeat]
                                    print 'bins'
                                    print bins
                                    print
                                    raise Exception('')
            
            #### two dimensional
            for l in indxpopl:
                for a, strgfrst in enumerate(liststrgfeatcorr[l]):
                    for b, strgseco in enumerate(liststrgfeatcorr[l]):
                        dictelem[l]['hist' + strgfrst + strgseco] = zeros((gdat.numbbinsplot, gdat.numbbinsplot))
                        if a < b:
                            binsfrst = getattr(gdat, 'bins' + strgfrst)
                            binsseco = getattr(gdat, 'bins' + strgseco)
                            if len(dictelem[l][strgfrst]) > 0 and len(dictelem[l][strgseco]) > 0:
                                dictelem[l]['hist' + strgfrst + strgseco] = histogram2d(dictelem[l][strgfrst][listindxfittelemfilt[0][l]], \
                                                                                             dictelem[l][strgseco][listindxfittelemfilt[0][l]], [binsfrst, binsseco])[0]
                            setattr(gdatobjt, strg + 'hist' + strgfrst + strgseco + 'pop%d' % l, dictelem[l]['hist' + strgfrst + strgseco])
            
            ### priors on element parameters and features
            for l in indxpopl:
                dictelem[l]['hist' + strgfeat + 'prio'] = empty(gdat.numbbinsplotprio)
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
                        dictelem[l]['hist' + strgfeat + 'prio'] = meanelem[l] * pdfn * deltprio * delt[0] / deltprio[0]
                    setattr(gdatobjt, strg + 'hist' + strgfeat + 'pop%dprio' % l, dictelem[l]['hist' + strgfeat + 'prio'])
                    if strg == 'true':
                        setattr(gdatobjt, 'refrhist' + strgfeat + 'pop%dprio' % l, dictelem[l]['hist' + strgfeat + 'prio'])
        
        # more derived parameters
        if gdat.elemtype == 'lens':
            if numbtrap > 0:
                if strgmodl != 'true': 
                    if gdat.priofactdoff != 0. or (gdat.fittnamefixp[gdat.fittindxfixpmeanelem] == 'logt').any():
                        # temp
                        indx = where(sum(gdat.refrhistmcut, axis=0) > 0)[0]
                        histmcutcorr = empty((numbpopl, gdat.numbbinsplot))
                        for l in indxpopl:
                            histmcutcorr[l, indx] = gdat.refrhistmcutpars[l, indx] * dictelem[l]['histmcut'][indx] / gdat.refrhistmcut[l, indx]
                        setattr(gdatobjt, strg + 'histmcutcorr', histmcutcorr)
            
                ## total truncated mass of the subhalo as a cross check
                # temp -- generalize
                if gdat.variasca:
                    asca = dictelem[0]['asca']
                else:
                    asca = gdat.ascaglob
                if gdat.variacut:
                    acut = dictelem[0]['acut']
                else:
                    acut = gdat.acutglob
                factmcutfromdefs = retr_factmcutfromdefs(gdat, gdat.adissour, gdat.adishost, gdat.adishostsour, asca, acut) 
                masssubh = array([sum(factmcutfromdefs * dictelem[0]['defs'])])

            ## calculate the host mass, subhalo mass and its fraction as a function of host halo-centric angle and interpolate at the Einstein radius of the host
            angl = sqrt((gdat.meanlgalcart - lgalhost)**2 + (gdat.meanbgalcart - bgalhost)**2)
            
            if lensmodltype != 'none':
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
                            dicttemp[name][k] = sum(convtemp[indxpixl]) * gdat.mdencrit * gdat.apix * gdat.adishost**2# / 2. / pi * gdat.deltanglhalf[k]
                        if name.endswith('intg'):
                            indxpixl = where(angl < gdat.meananglhalf[k])
                            dicttemp[name][k] = sum(convtemp[indxpixl]) * gdat.mdencrit * gdat.apix * gdat.adishost**2# * indxpixl[0].size
                        
                        if name[:4] == 'frac':
                            if dicttemp['masshost' + name[8:]][k] != 0.:
                                dicttemp['fracsubh' + name[8:]][k] = dicttemp['masssubh' + name[8:]][k] / dicttemp['masshost' + name[8:]][k]
                    setattr(gdatobjt, strg + name, dicttemp[name])
                
                # interpolate the host, subhalo masses and subhalo mass fraction at the Einstein radius and save it as a scalar variable
                for name in listnamevarbmassvect:
                    dicttemp[name + 'bein'] = interp(beinhost, gdat.meananglhalf, dicttemp[name])
                    setattr(gdatobjt, strg + name + 'bein', array([dicttemp[name + 'bein']]))
            
        if numbtrap > 0:
            ## copy element features to the global object
            feat = [[] for l in indxpopl]
            for l in indxpopl:
                feat[l] = dict()
                for strgfeat in liststrgfeat[l]:
                    if strg == 'true':
                        shap = list(ones(dictelem[l][strgfeat].ndim, dtype=int))
                        feat[l][strgfeat] = tile(dictelem[l][strgfeat], [3] + shap)
                    if strg == 'this':
                        feat[l][strgfeat] = dictelem[l][strgfeat]
                    setattr(gdatobjt, strg + 'hist' + strgfeat + 'pop%d' % l, dictelem[l]['hist' + strgfeat])
                    if strg == 'true':
                        setattr(gdatobjt, 'refrhist' + strgfeat + 'pop%d' % l, dictelem[l]['hist' + strgfeat])
            for strgfeat in liststrgfeattotl:
                feattemp = [[] for l in indxpopl]
                for l in indxpopl:
                    if strgfeat in feat[l]:
                        feattemp[l] = feat[l][strgfeat]
                    else:
                        feattemp[l] = array([])
                if strg == 'true':
                    liststrgtemp = ['true', 'refr']
                else:
                    liststrgtemp = [strg]
                for strgtemp in liststrgtemp:
                    setattr(gdatobjt, strgtemp + strgfeat, feattemp)
                
        ### Exculusive comparison with the true state
        if strg != 'true' and gdat.datatype == 'mock':
            
            if gdat.elemtype == 'lens' and gdat.datatype == 'mock':
                deflsingresi = deflsing - gdat.truedeflsing
                deflresi = defl - gdat.truedefl
                deflmgtdresi = sqrt(sum(deflresi**2, axis=2))
                deflresiperc = 100. * deflmgtdresi / gdat.truedeflmgtd[:, :, None]
                deflsingresimgtd = sqrt(sum(deflsingresi**2, axis=2))
                deflsingresiperc = 100. * deflsingresimgtd / gdat.truedeflsingmgtd
                setattr(gdatobjt, strg + 'deflsingresi', deflsingresi)
                setattr(gdatobjt, strg + 'deflresi', deflresi)
                setattr(gdatobjt, strg + 'deflmgtdresi', deflmgtdresi)
                if numbtrap > 0:
                    convelemresi = convelem - gdat.trueconvelem
                    convelemresiperc = 100. * convelemresi / gdat.trueconvelem
                    setattr(gdatobjt, strg + 'convelemresi', convelemresi)
                    setattr(gdatobjt, strg + 'convelemresiperc', convelemresiperc)
                magnresi = magn - gdat.truemagn
                magnresiperc = 100. * magnresi / gdat.truemagn
                setattr(gdatobjt, strg + 'magnresi', magnresi)
                setattr(gdatobjt, strg + 'magnresiperc', magnresiperc)
    
        if numbtrap > 0:
            # correlate the catalog sample with the reference catalog
            if gdat.refrinfo and gdat.refrlgal != None and gdat.refrbgal != None and not (strg == 'true' and gdat.datatype == 'mock'):
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        if gdat.refrnumbelem[q] > 0:
                            cmpl = array([float(len(indxelemrefrasschits[q][l])) / gdat.refrnumbelem[q]])
                            if gdat.diagmode:
                                if cmpl > 1. or cmpl < 0.:
                                    raise Exception('')
                        else:
                            cmpl = -1.
                        setattr(gdatobjt, strg + 'cmplref%dpop%d' % (q, l), cmpl)
                        if numbelem[l] > 0:
                            fdis = array([float(indxelemfittasscfals[q][l].size) / numbelem[l]])
                            if gdat.diagmode:
                                if fdis > 1. or fdis < 0.:
                                    raise Exception('')
                        else:
                            fdis = -1.
                        setattr(gdatobjt, strg + 'fdisref%dpop%d' % (q, l), fdis)
                
                # collect the associated fitting element feature for each reference element
                featrefrassc = [[[] for l in indxpopl] for q in gdat.indxrefr]
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        featrefrassc[q][l] = dict()
                        for strgfeat in gdat.listnamefeatrefr[q]:
                            if strgfeat.endswith('pars') or strgfeat.endswith('nrel') or not strgfeat in liststrgfeat[l]:
                                continue
                            if isinstance(dictelem[l][strgfeat], ndarray) and dictelem[l][strgfeat].ndim > 1:
                                continue
                            featrefrassc[q][l][strgfeat] = zeros(gdat.refrnumbelem[q]) 
                            
                            if len(indxelemrefrasschits[q][l]) > 0 and len(dictelem[l][strgfeat]) > 0:
                                featrefrassc[q][l][strgfeat][indxelemrefrasschits[q][l]] = dictelem[l][strgfeat][indxelemfittasschits[q][l]]
                            setattr(gdatobjt, strg + strgfeat + 'asscref%dpop%d' % (q, l), featrefrassc[q][l][strgfeat])
                
                # completeness
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        for strgfeat in liststrgfeatodim[l]:
                            cmplfeat = zeros(gdat.numbbinsplot)
                            errrcmplfeat = zeros((2, gdat.numbbinsplot))
                            refrfeat = getattr(gdat, 'refr' + strgfeat)
                            if refrfeat == None or not strgfeat in liststrgfeatodim[q]:
                                continue
                            refrhistfeat = getattr(gdat, 'refrhist' + strgfeat)
                            #indx = where(isfinite(featrefrassc[q][l][strgfeat]))[0]
                            bins = getattr(gdat, 'bins' + strgfeat)
                            if len(indxelemrefrasschits[q][l]) > 0:
                                histfeatrefrassc = histogram(refrfeat[q][0, indxelemrefrasschits[q][l]], bins=bins)[0]
                                if refrhistfeat != None:
                                    cmplfeat = histfeatrefrassc / refrhistfeat[q, :]
                                    errrcmplfeat[:, :] = (cmplfeat / sqrt(maximum(ones(gdat.numbbinsplot), refrhistfeat[q, :])))[None, :]
                                    if gdat.diagmode:
                                        if where((cmplfeat > 1.) | (cmplfeat < 0.))[0].size > 0:
                                            raise Exception('')
                            setattr(gdatobjt, strg + 'cmpl' + strgfeat + 'ref%dpop%d' % (q, l), cmplfeat)
                            setattr(gdatobjt, strg + 'errrcmpl' + strgfeat + 'ref%dpop%d' % (q, l), errrcmplfeat)
                   
                # false discovery rate
                for q in gdat.indxrefr:
                    for l in gdat.fittindxpopl:
                        for strgfeat in liststrgfeatodim[l]:
                            fdisfeat = zeros(gdat.numbbinsplot)
                            errrfdisfeat = zeros((2, gdat.numbbinsplot))
                            if not strgfeat in liststrgfeatodim[l] or strgfeat in gdat.listnamefeatrefronly[q][l]:
                                continue
                            bins = getattr(gdat, 'bins' + strgfeat)
                            if len(indxelemfittasscfals[q][l]) > 0 and len(dictelem[l][strgfeat]) > 0:
                                histfeatfals = histogram(dictelem[l][strgfeat][indxelemfittasscfals[q][l]], bins=bins)[0]
                                fitthistfeat = getattr(gdatobjt, strg + 'hist' + strgfeat + 'pop%d' % l)
                                fdisfeat = histfeatfals / fitthistfeat
                                errrfdisfeat[:, :] = (fdisfeat / sqrt(maximum(ones(gdat.numbbinsplot), fitthistfeat)))[None, :]
                                if gdat.diagmode:
                                    if where((fdisfeat > 1.) | (fdisfeat < 0.))[0].size > 0:
                                        raise Exception('')
                            setattr(gdatobjt, strg + 'fdis' + strgfeat + 'ref%dpop%d' % (q, l), fdisfeat)
                            setattr(gdatobjt, strg + 'errrfdis' + strgfeat + 'ref%dpop%d' % (q, l), errrfdisfeat)

            # temp
            if strgmodl == 'true' and gdat.verbtype > 0:
                for l in indxpopl:
                    for strgfeat in liststrgfeat[l]:
                        minm = getattr(gdat, 'minm' + strgfeat)
                        maxm = getattr(gdat, 'maxm' + strgfeat)
                        if where(minm > dictelem[l][strgfeat])[0].size > 0 or where(maxm < dictelem[l][strgfeat])[0].size > 0:
                            print 'Warning: element feature outside the plot limits.'
                            print 'Feature: '
                            print strgfeat
                            print 'Plot minmimum'
                            print minm
                            print 'Plot maxmimum'
                            print maxm
                            print 'Feature minimum'
                            print amin(dictelem[l][strgfeat])
                            print 'Feature maximum'
                            print amax(dictelem[l][strgfeat])
                            print
                            if strgfeat == gdat.namecompampl:
                                raise Exception('')


def proc_cntpdata(gdat):

    # exclude voxels with vanishing exposure
    ## data counts
    if gdat.datatype == 'inpt':
        gdat.cntpdata = gdat.cntpdata[gdat.indxcuberofi]
    
    ## spatial average
    gdat.sbrtdatamean = retr_spatmean(gdat, gdat.cntpdata, boolcntp=True)
    
    gdat.sbrtdatabrod = sum(sum(gdat.cntpdata, axis=2), axis=0)

    # obtain cartesian versions of the maps
    if gdat.pixltype == 'cart':
        ## data counts
        cntpdatacarttemp = zeros((gdat.numbener, gdat.numbpixlfull, gdat.numbevtt))
        cntpdatacarttemp[:, gdat.indxpixlrofi, :] = gdat.cntpdata
        gdat.cntpdatacart = cntpdatacarttemp.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
   

def retr_info(pdfnpost, pdfnprio):
    
    info = pdfnpost * log(pdfnpost / pdfnprio)

    return info


def retr_llik_depr(gdat, cntpmodl, indxpixleval):
    
    if gdat.liketype == 'pois':
    	llik = gdat.cntpdata[:, indxpixleval, :] * log(cntpmodl) - cntpmodl
    if gdat.liketype == 'gaus':
        llik = -0.5 * (gdat.cntpdata[:, indxpixleval, :] - cntpmodl[:, indxpixleval, :])**2 / gdat.cntpdata[:, indxpixleval, :]
    
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
        sbrtsers = spec[:, 0, None, None] * sp.interpolate.interpn((gdat.binslgalsers, gdat.binshalfsers, gdat.binsindxsers), gdat.sersprof, inpt)[None, :, None]
    
    # evaluate directly de Vaucouleurs
    if gdat.serstype == 'vauc':
        sbrtsers = swapaxes(spec[:, 0, None, None] * retr_sbrtsersnorm(angl, size)[None, :, None], 1, 2)
        sbrtsers = spec[:, 0, None, None] * retr_sbrtsersnorm(angl, size)[None, :, None]
    
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


def prep_gdatmodi(gdat, gdatmodi, gdatobjt, strg):
    
    if gdat.verbtype > 1:
        print 'prep_gdatmodi()'

    gdatmodi.thischro = zeros(gdat.numbchro)
    gdatmodi.thistmprfactstdv = 1.
    for namevarbstat in gdat.listnamevarbstat:
        temp = copytdgu(getattr(gdatobjt, strg + namevarbstat))
        setattr(gdatmodi, 'this' + namevarbstat, temp)
    

def make_anim(strgsrch, full=False):

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
    

