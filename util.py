# common imports
from __init__ import *

class globdatastrt(object):
    
    def __init__(self):
        pass


def retr_psfnwdth(globdata, psfn, frac):
    
    fwhm = zeros((globdata.numbener, globdata.numbevtt))
    for i in globdata.indxener:
        for m in globdata.indxevtt:
            indxangldispgood = argsort(psfn[i, :, m])
            intpfwhm = max(frac * amax(psfn[i, :, m]), amin(psfn[i, :, m]))
            if intpfwhm > amin(psfn[i, indxangldispgood, m]) and intpfwhm < amax(psfn[i, indxangldispgood, m]):
                fwhm[i, m] = interp1d(psfn[i, indxangldispgood, m], globdata.angldisp[indxangldispgood])(intpfwhm)
    return fwhm


def retr_spec(globdata, flux, sind):

    if isscalar(flux):
        flux = array([flux])

    if isscalar(sind):
        sind = array([sind])
        
    spec = flux[None, :] * (globdata.meanener[:, None] / globdata.enerfdfn)**(-sind[None, :])
    
    return spec


def retr_indx(globdata, indxpntsfull):    

    indxsamplgal = []
    indxsampbgal = []
    indxsampspec = []
    indxsampsind = []
    indxsampcomp = []
    for l in globdata.indxpopl:
        indxsamplgaltemp = globdata.indxsampcompinit + globdata.maxmnumbcomp * l +             array(indxpntsfull[l], dtype=int) * globdata.numbcomp
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], globdata.numbener, 0) +repeat(arange(globdata.numbener), len(indxpntsfull[l])).reshape(globdata.numbener, -1))
        if globdata.colrprio:
            indxsampsind.append(indxsamplgaltemp + 2 + globdata.numbener)
        indxsampcomp.append(repeat(indxsamplgaltemp, globdata.numbcomp) + tile(arange(globdata.numbcomp, dtype=int), len(indxpntsfull[l])))

    return indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp


def retr_pntsflux(globdata, lgal, bgal, spec, psfipara, psfntype):
    
    numbpnts = lgal.size
    
    # calculate the distance to all pixels from each point source
    dist = empty((globdata.numbpixl, numbpnts))
    for k in range(numbpnts):
        dist[:, k] = retr_dist(globdata, lgal[k], bgal[k], globdata.lgalgrid, globdata.bgalgrid)

    # evaluate the PSF
    pntsflux = empty((numbpnts, globdata.numbener, globdata.numbpixl, globdata.numbevtt))
    for k in range(numbpnts):
        psfn = retr_psfn(globdata, psfipara, globdata.indxener, dist[:, k], psfntype)
        pntsflux[k, :, :, :] = spec[:, k, None, None] * psfn

    # sum contributions from all PS
    pntsfluxtemp = sum(pntsflux, 0) 

    return pntsfluxtemp


def retr_rofi_flux(globdata, normback, pntsflux, tempindx):

    modlflux = pntsflux[tempindx]
    for c in globdata.indxback:
        modlflux += normback[c, :, None, None] * globdata.backflux[c][tempindx]        
    
    return modlflux


def cdfn_spec_brok(globdata, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):

    norm = 1. / ((1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) +                  ((globdata.maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))
    if flux <= fluxbrek:
        fluxunit = norm / (1. - fdfnsloplowr) *             ((flux / fluxbrek)**(1. - fdfnsloplowr) - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))
    else:
        fluxunit = norm / (1. - fdfnsloplowr) * (1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) -             norm / (1. - fdfnslopuppr) * (1. - (flux / fluxbrek)**(1. - fdfnslopuppr))
       
    return fluxunit


def pdfn_spec_brok(globdata, flux, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):

    norm = 1. / ((1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) +                  ((globdata.maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr) - 1.) / (1. - fdfnslopuppr))
    if flux <= fluxbrek:
        pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnsloplowr)
    else:
        pdfnflux = norm * (flux / fluxbrek)**(1. - fdfnslopuppr)
        
    return pdfnflux


def icdf_spec_brok(globdata, fluxunit, fdfnsloplowr, fdfnslopuppr, fluxbrek, i):
    
    norm = 1. / ((1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) / (1. - fdfnsloplowr) -                  (1. - (globdata.maxmspec[i] / fluxbrek)**(1. - fdfnslopuppr)) / (1. - fdfnslopuppr))
    fluxunitbrek = norm / (1. - fdfnsloplowr) * (1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))
    if fluxunit < fluxunitbrek:
        flux = fluxbrek * (fluxunit * (1. - fdfnsloplowr) / norm +                        (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr))**(1. / (1. - fdfnsloplowr))
    else:
        flux = fluxbrek * (1. - (norm / (1. - fdfnsloplowr) * (1. - (globdata.minmspec[i] / fluxbrek)**(1. - fdfnsloplowr)) -        fluxunit) * (1. - fdfnslopuppr) / norm)**(1. / (1. - fdfnslopuppr))

    return flux


def cdfn_spec(globdata, flux, fdfnslop, minmspectemp, maxmspectemp):
        
    fluxunit = (flux**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop)) / (maxmspectemp**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop))
        
    return fluxunit


def icdf_spec(globdata, fluxunit, fdfnslop, minmspectemp, maxmspectemp):
    
    flux = (fluxunit * (maxmspectemp**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop)) + minmspectemp**(1. - fdfnslop))**(1. / (1. - fdfnslop))
    
    return flux


def pdfn_spec(globdata, flux, fdfnslop, minmspectemp, maxmspectemp):
  
    pdfnflux = (1. - fdfnslop) / (maxmspectemp**(1. - fdfnslop) - minmspectemp**(1. - fdfnslop)) * flux**(-fdfnslop)
          
    return pdfnflux


def icdf_self(paraunit, minmpara, factpara):
    para = factpara * paraunit + minmpara
    return para


def cdfn_self(para, minmpara, factpara):
    paraunit = (para - minmpara) / factpara
    return paraunit


def icdf_logt(paraunit, minmpara, factpara):
    para = minmpara * exp(paraunit * factpara)
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


def icdf_psfipara(globdata, psfiparaunit, thisindxpsfipara):

    minmpsfipara = globdata.minmpsfipara[thisindxpsfipara]
    factpsfipara = globdata.factpsfipara[thisindxpsfipara]
    scalpsfipara = globdata.scalpsfipara[thisindxpsfipara]
        
    if scalpsfipara == 'self':
        psfipara = icdf_self(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfipara = icdf_logt(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfipara = icdf_atan(psfiparaunit, minmpsfipara, factpsfipara)

    return psfipara


def cdfn_psfipara(globdata, psfipara, thisindxpsfipara):
    
    minmpsfipara = globdata.minmpsfipara[thisindxpsfipara]
    factpsfipara = globdata.factpsfipara[thisindxpsfipara]
    scalpsfipara = globdata.scalpsfipara[thisindxpsfipara]

    if scalpsfipara == 'self':
        psfiparaunit = cdfn_self(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfiparaunit = cdfn_logt(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfiparaunit = cdfn_atan(psfipara, minmpsfipara, factpsfipara)
    
    return psfiparaunit


def retr_thisindxprop(globdata, samp):
    
    # choose the population to be modified
    globdata.indxpoplmodi = choice(globdata.indxpopl)
    
    numbpnts = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]]
    if numbpnts == globdata.maxmnumbpnts[globdata.indxpoplmodi]:
        globdata.thisindxprop = choice(globdata.indxprop, p=globdata.probpropmaxm)
    elif numbpnts == globdata.minmnumbpnts:
        globdata.thisindxprop = choice(globdata.indxprop, p=globdata.probpropminm)
    else:
        globdata.thisindxprop = choice(globdata.indxprop, p=globdata.probprop)
        

def retr_indxpixl(globdata, bgal, lgal):

    if globdata.pixltype == 'heal':
        indxpixl = globdata.pixlcnvt[ang2pix(globdata.numbsideheal, deg2rad(90. - bgal), deg2rad(lgal))]
    else:
        
        indxlgcr = floor(globdata.numbsidecart * (lgal - globdata.minmlgal) / 2. / globdata.maxmgang).astype(int)
        indxbgcr = floor(globdata.numbsidecart * (bgal - globdata.minmbgal) / 2. / globdata.maxmgang).astype(int)

        if isscalar(indxlgcr):
            if indxlgcr < 0:
                indxlgcr = 0
            if indxlgcr >= globdata.numbsidecart:
                indxlgcr = globdata.numbsidecart - 1
        else:
            indxlgcr[where(indxlgcr < 0)] = 0
            indxlgcr[where(indxlgcr >= globdata.numbsidecart)] = globdata.numbsidecart - 1
            
        if isscalar(indxbgcr):
            if indxbgcr < 0:
                indxbgcr = 0
            if indxbgcr >= globdata.numbsidecart:
                indxbgcr = globdata.numbsidecart - 1
        else:
            indxbgcr[where(indxbgcr < 0)] = 0
            indxbgcr[where(indxbgcr >= globdata.numbsidecart)] = globdata.numbsidecart - 1
            
        indxpixl = indxbgcr * globdata.numbsidecart + indxlgcr

    return indxpixl


def retr_llik(globdata, init=False):

    if init:

        globdata.thisllik = globdata.datacnts * log(globdata.thismodlcnts) - globdata.thismodlcnts
        
    elif globdata.thisindxprop >= globdata.indxproppsfipara:

        # determine pixels over which to evaluate the log-likelihood
        if globdata.thisindxprop == globdata.indxpropnormback:
            indxpixlmodi = globdata.indxpixl
            
        if globdata.thisindxprop == globdata.indxproppsfipara or globdata.thisindxprop >= globdata.indxpropbrth:
            
            if globdata.thisindxprop == globdata.indxproppsfipara:
                numbpnts = int(sum(globdata.thissampvarb[globdata.indxsampnumbpnts]))
                lgal = globdata.thissampvarb[concatenate(globdata.thisindxsamplgal)]
                bgal = globdata.thissampvarb[concatenate(globdata.thisindxsampbgal)]
                spec = globdata.thissampvarb[concatenate(globdata.thisindxsampspec)[globdata.indxenermodi, :]]
            else:
                numbpnts = globdata.modilgal.shape[0]
                lgal = globdata.modilgal
                bgal = globdata.modibgal
                spec = globdata.modispec
                
            thisindxpixlprox = []
            for k in range(numbpnts):
                indxspecproxtemp = argmin(spec[0, k] - globdata.meanspecprox)
                indxpixltemp = retr_indxpixl(globdata, bgal[k], lgal[k])
                thisindxpixlprox.append(globdata.indxpixlprox[indxspecproxtemp][indxpixltemp])
            indxpixlmodi = unique(concatenate(thisindxpixlprox))

        # construct the mesh grid for likelihood evaluation
        if globdata.thisindxprop >= globdata.indxproppsfipara:
            globdata.indxcubemodi = meshgrid(globdata.indxenermodi, indxpixlmodi, globdata.indxevtt, indexing='ij')

        # update the point source flux map
        if globdata.thisindxprop == globdata.indxproppsfipara or globdata.thisindxprop >= globdata.indxpropbrth:

            if globdata.thisindxprop == globdata.indxproppsfipara:
                globdata.nextpntsflux[globdata.indxcubemodi] = 0.
            else:
                globdata.nextpntsflux[globdata.indxcubemodi] = globdata.thispntsflux[globdata.indxcubemodi]
                
            for k in range(numbpnts):
                
                # calculate the distance to the pixels to be updated
                dist = retr_dist(globdata, lgal[k], bgal[k], globdata.lgalgrid[thisindxpixlprox[k]], globdata.bgalgrid[thisindxpixlprox[k]])

                # evaluate the PSF over the set of data cubes to be updated
                temppsfipara = copy(globdata.thissampvarb[globdata.indxsamppsfipara])
                if globdata.thisindxprop == globdata.indxproppsfipara:
                    temppsfipara[globdata.indxpsfiparamodi] = globdata.nextsampvarb[globdata.indxsampmodi]
                psfn = retr_psfn(globdata, temppsfipara, globdata.indxenermodi, dist, globdata.psfntype)
      
                # update the data cubes
                for i in range(globdata.indxenermodi.size):
                    globdata.nextpntsflux[globdata.indxenermodi[i], thisindxpixlprox[k], :] += spec[i, k] * psfn[i, :, :]
            
        if globdata.thisindxprop == globdata.indxpropnormback:
            
            normbacktemp = empty((globdata.numbback, 1))
            for c in globdata.indxback:
                if c == globdata.indxbackmodi:
                    normbacktemp[c, 0] = globdata.nextsampvarb[globdata.indxsampnormback[c, globdata.indxenermodi]]
                else:
                    normbacktemp[c, 0] = globdata.thissampvarb[globdata.indxsampnormback[c, globdata.indxenermodi]]
                
            globdata.nextmodlflux[globdata.indxcubemodi] = retr_rofi_flux(globdata, normbacktemp, globdata.thispntsflux, globdata.indxcubemodi)

        if globdata.thisindxprop == globdata.indxproppsfipara or globdata.thisindxprop >= globdata.indxpropbrth:
            normback = globdata.thissampvarb[globdata.indxsampnormback[meshgrid(globdata.indxback, globdata.indxenermodi, indexing='ij')]]
            globdata.nextmodlflux[globdata.indxcubemodi] = retr_rofi_flux(globdata, normback, globdata.nextpntsflux, globdata.indxcubemodi)

        globdata.nextmodlcnts[globdata.indxcubemodi] = globdata.nextmodlflux[globdata.indxcubemodi] * globdata.expo[globdata.indxcubemodi] * \
            globdata.apix * globdata.diffener[globdata.indxenermodi, None, None] # [1]
        globdata.nextllik[globdata.indxcubemodi] = globdata.datacnts[globdata.indxcubemodi] * log(globdata.nextmodlcnts[globdata.indxcubemodi]) - \
            globdata.nextmodlcnts[globdata.indxcubemodi]
            
        if not isfinite(globdata.nextllik[globdata.indxcubemodi]).any():
            warnings.warn('Log-likelihood went NAN!')
            
        globdata.deltllik = sum(globdata.nextllik[globdata.indxcubemodi] - globdata.thisllik[globdata.indxcubemodi])
    else:
        globdata.deltllik = 0.
        
    
def retr_fdfn(globdata, fdfnnorm, fdfnslop, i):
               
    fluxhistmodl = fdfnnorm / globdata.diffspec[i, globdata.indxspecpivt] *         globdata.diffspec[i, :] * (globdata.meanspec[i, :] / globdata.fluxpivt[i])**(-fdfnslop)
          
    return fluxhistmodl


def retr_lpri(globdata, init=False):
        
    if init:
        globdata.thislpri = zeros((globdata.numbpopl, globdata.numbener))
        
        for i in globdata.indxenerfdfn:
            for l in globdata.indxpopl:
                fdfnnorm = globdata.thissampvarb[globdata.indxsampfdfnnorm[l]]
                fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[l, i]]
                fluxhistmodl = retr_fdfn(globdata, fdfnnorm, fdfnslop, i)
                spec = globdata.thissampvarb[globdata.thisindxsampspec[l][i, :]]
                fluxhist = histogram(spec, globdata.binsspec[i, :])[0]
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                globdata.thislpri[l, i] = sum(lprbpois)            
      
        globdata.nextlpri = copy(globdata.thislpri)
                
    else:
        globdata.nextlpri = copy(globdata.thislpri)
        if globdata.thisindxprop == globdata.indxpropfdfnnorm or globdata.thisindxprop == globdata.indxpropfdfnslop             or globdata.thisindxprop >= globdata.indxpropbrth and globdata.thisindxprop <= globdata.indxpropmerg:
              
            if globdata.colrprio:
                indxenertemp = globdata.indxenerfdfn
            else:
                indxenertemp = globdata.indxenermodi
            for i in indxenertemp:
                if globdata.thisindxprop == globdata.indxpropfdfnnorm:
                    fdfnnorm = globdata.nextsampvarb[globdata.indxsampfdfnnorm[globdata.indxpoplmodi]]
                    fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]]
                elif globdata.thisindxprop == globdata.indxpropfdfnslop:
                    fdfnnorm = globdata.thissampvarb[globdata.indxsampfdfnnorm[globdata.indxpoplmodi]]
                    fdfnslop = globdata.nextsampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]]
                else:
                    fdfnnorm = globdata.thissampvarb[globdata.indxsampfdfnnorm[globdata.indxpoplmodi]]
                    fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]]
                fluxhistmodl = retr_fdfn(globdata, fdfnnorm, fdfnslop, i)
    
                fluxhist = histogram(globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][i, :]], globdata.binsspec[i, :])[0] 
                if globdata.thisindxprop == globdata.indxpropbrth:
                    fluxhist += histogram(globdata.modispec[i, 0], globdata.binsspec[i, :])[0]
                elif globdata.thisindxprop == globdata.indxpropdeth:
                    fluxhist -= histogram(-globdata.modispec[i, 0], globdata.binsspec[i, :])[0]
                elif globdata.thisindxprop == globdata.indxpropsplt:
                    fluxhist -= histogram(-globdata.modispec[i, 0], globdata.binsspec[i, :])[0]
                    fluxhist += histogram(globdata.modispec[i, 1:3], globdata.binsspec[i, :])[0]
                elif globdata.thisindxprop == globdata.indxpropmerg:
                    fluxhist -= histogram(-globdata.modispec[i, 0:2], globdata.binsspec[i, :])[0]
                    fluxhist += histogram(globdata.modispec[i, 2], globdata.binsspec[i, :])[0]
                
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                globdata.nextlpri[globdata.indxpoplmodi, i] = sum(lprbpois)

            globdata.deltlpri = sum(globdata.nextlpri[globdata.indxpoplmodi, globdata.indxenermodi] - globdata.thislpri[globdata.indxpoplmodi, globdata.indxenermodi])
        else:
            globdata.deltlpri = 0.

        
def pars_samp(globdata, indxpntsfull, samp):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(globdata, indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[globdata.indxsampnumbpnts] = samp[globdata.indxsampnumbpnts]
    sampvarb[globdata.indxsampfdfnnorm] = icdf_logt(samp[globdata.indxsampfdfnnorm], globdata.minmfdfnnorm, globdata.factfdfnnorm)
    sampvarb[globdata.indxsampfdfnslop] = icdf_atan(samp[globdata.indxsampfdfnslop], globdata.minmfdfnslop, globdata.factfdfnslop)
    for c in globdata.indxback:
        sampvarb[globdata.indxsampnormback[c, :]] = icdf_logt(samp[globdata.indxsampnormback[c, :]], globdata.minmnormback[c], globdata.factnormback[c])
    for k in globdata.indxpsfipara:
        sampvarb[globdata.indxsamppsfipara[k]] = icdf_psfipara(globdata, samp[globdata.indxsamppsfipara[k]], k)
    cnts = []
    listspectemp = []
    indxpixlpnts = []
    for l in globdata.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg) 
        for i in globdata.indxenerfdfn:
            sampvarb[indxsampspec[l][i, :]] = icdf_spec(globdata, samp[indxsampspec[l][i, :]], sampvarb[globdata.indxsampfdfnslop[l, i]], globdata.minmspec[i], globdata.maxmspec[i])
        if globdata.colrprio:
            sampvarb[indxsampsind[l]] = icdf_atan(samp[indxsampsind[l]], globdata.minmsind, globdata.factsind)
            sampvarb[indxsampspec[l]] = retr_spec(globdata, sampvarb[indxsampspec[l][globdata.indxenerfdfn, :]], sampvarb[indxsampsind[l]])
            
        indxpixlpntstemp = retr_indxpixl(globdata, sampvarb[indxsampbgal[l]], sampvarb[indxsamplgal[l]])
    
        cntstemp = sampvarb[indxsampspec[l]][:, :, None] * globdata.expo[:, indxpixlpntstemp, :] * globdata.diffener[:, None, None]
            
        cnts.append(cntstemp)
        indxpixlpnts.append(indxpixlpntstemp)
        listspectemp.append(sampvarb[indxsampspec[l]])
        
    pntsflux = retr_pntsflux(globdata,    sampvarb[concatenate(indxsamplgal)],    sampvarb[concatenate(indxsampbgal)],    concatenate(listspectemp, axis=1),    sampvarb[globdata.indxsamppsfipara], globdata.psfntype)
    
    totlflux = retr_rofi_flux(globdata, sampvarb[globdata.indxsampnormback], pntsflux, globdata.indxcubefull)
    totlcnts = totlflux * globdata.expo * globdata.apix * globdata.diffener[:, None, None] # [1]
    
    return sampvarb, indxpixlpnts, cnts, pntsflux, totlflux, totlcnts
    
    
def retr_mrkrsize(globdata, spec, indxenertemp):

    mrkrsize = (spec - globdata.minmspec[indxenertemp]) / (globdata.maxmspec[indxenertemp] - globdata.minmspec[indxenertemp]) * \
		(globdata.maxmmrkrsize - globdata.minmmrkrsize) + globdata.minmmrkrsize
        
    return mrkrsize


def retr_scalfromangl(globdata, thisangl, i, m):
    scalangl = 2. * arcsin(.5 * sqrt(2. - 2 * cos(thisangl))) / globdata.fermscalfact[i, m]
    return scalangl

def retr_anglfromscal(globdata, scalangl, i, m):
    thisangl = arccos(1. - 0.5 * (2. * sin(scalangl * globdata.fermscalfact[i, m] / 2.))**2)
    return thisangl


def retr_fermpsfn(globdata):
   
    name = os.environ["PCAT_DATA_PATH"] + '/irf/psf_P8R2_SOURCE_V6_PSF.fits'
    irfn = pf.getdata(name, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = sqrt(minmener * maxmener)

    parastrg = ['ntail', 'score', 'gcore', 'stail', 'gtail']

    globdata.numbfermformpara = 5
    globdata.numbfermscalpara = 3
    
    fermscal = zeros((globdata.numbevtt, globdata.numbfermscalpara))
    fermform = zeros((globdata.numbener, globdata.numbevtt, globdata.numbfermformpara))
    globdata.fermpsfipara = zeros((globdata.numbener * globdata.numbfermformpara * globdata.numbevtt))
    
    for m in globdata.indxevtt:
        fermscal[m, :] = pf.getdata(name, 2 + 3 * globdata.indxevttincl[m])['PSFSCALE']
        irfn = pf.getdata(name, 1 + 3 * globdata.indxevttincl[m])
        for k in range(5):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(globdata.meanener)

    globdata.fermscalfact = sqrt((fermscal[None, :, 0] * (10. * globdata.meanener[:, None])**fermscal[None, :, 2])**2 +                         fermscal[None, :, 1]**2)
    
    # convert N_tail to f_core
    for m in globdata.indxevtt:
        for i in globdata.indxener:
            fermform[i, m, 1] = retr_anglfromscal(globdata, fermform[i, m, 1], i, m) # [rad]
            fermform[i, m, 3] = retr_anglfromscal(globdata, fermform[i, m, 3], i, m) # [rad]
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)
    
    # store the fermi PSF parameters
    for m in globdata.indxevtt:
        for k in range(globdata.numbfermformpara):
            indxfermpsfiparatemp = m * globdata.numbfermformpara * globdata.numbener + globdata.indxener * globdata.numbfermformpara + k
            globdata.fermpsfipara[indxfermpsfiparatemp] = fermform[:, m, k]
        
    frac = fermform[:, :, 0]
    sigc = fermform[:, :, 1]
    gamc = fermform[:, :, 2]
    sigt = fermform[:, :, 3]
    gamt = fermform[:, :, 4]
    
    globdata.fermpsfn = retr_doubking(globdata.angldisp[None, :, None], frac[:, None, :], sigc[:, None, :], gamc[:, None, :],sigt[:, None, :], gamt[:, None, :])


def updt_samp(globdata):
    
    if globdata.thisindxprop == globdata.indxpropfdfnnorm:
        globdata.thissampvarb[globdata.indxsampfdfnnorm[globdata.indxpoplmodi]] = globdata.nextsampvarb[globdata.indxsampfdfnnorm[globdata.indxpoplmodi]]
        globdata.thislpri[globdata.indxpoplmodi, globdata.indxenermodi] = globdata.nextlpri[globdata.indxpoplmodi, globdata.indxenermodi]

    if globdata.thisindxprop == globdata.indxpropfdfnslop:
 
        # update the unit sample vector
        globdata.drmcsamp[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, :], -1] = \
            cdfn_spec(globdata, globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, :]], \
            globdata.nextsampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]], \
            globdata.minmspec[globdata.indxenermodi], globdata.maxmspec[globdata.indxenermodi])
        
        # update the sample vector
        globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]] = \
            globdata.nextsampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]]
            
        # update the prior register
        globdata.thislpri[globdata.indxpoplmodi, globdata.indxenermodi] = globdata.nextlpri[globdata.indxpoplmodi, globdata.indxenermodi]

    # likelihood updates
    if globdata.thisindxprop >= globdata.indxproppsfipara:
        globdata.thisllik[globdata.indxcubemodi] = globdata.nextllik[globdata.indxcubemodi]
        globdata.thismodlcnts[globdata.indxcubemodi] = globdata.nextmodlcnts[globdata.indxcubemodi]
        
    if globdata.thisindxprop == globdata.indxproppsfipara:
        globdata.thissampvarb[globdata.indxsampmodi] = globdata.nextsampvarb[globdata.indxsampmodi]
        
    if globdata.thisindxprop == globdata.indxpropnormback:
        globdata.thissampvarb[globdata.indxsampnormback[globdata.indxbackmodi, globdata.indxenermodi]] = \
            globdata.nextsampvarb[globdata.indxsampnormback[globdata.indxbackmodi, globdata.indxenermodi]]
        
    if globdata.thisindxprop >= globdata.indxpropbrth or globdata.thisindxprop == globdata.indxproppsfipara:
        globdata.thispntsflux[globdata.indxcubemodi] = globdata.nextpntsflux[globdata.indxcubemodi]
        
    # transdimensinal updates
    if globdata.thisindxprop >= globdata.indxpropbrth and globdata.thisindxprop <= globdata.indxpropmerg:
        globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]]
        globdata.thislpri[globdata.indxpoplmodi, globdata.indxenermodi] = globdata.nextlpri[globdata.indxpoplmodi, globdata.indxenermodi]
        
    # birth
    if globdata.thisindxprop == globdata.indxpropbrth:
        
        # update the PS index lists
        globdata.thisindxpntsfull[globdata.indxpoplmodi].append(globdata.thisindxpntsempt[globdata.indxpoplmodi][0])
        del globdata.thisindxpntsempt[globdata.indxpoplmodi][0]

        # update the components
        globdata.thissampvarb[globdata.indxsampmodi[0]] = globdata.modilgal
        globdata.thissampvarb[globdata.indxsampmodi[1]] = globdata.modibgal
        if globdata.colrprio:
            globdata.thissampvarb[globdata.indxsampmodi[2+globdata.indxener]] = globdata.modispec
            globdata.thissampvarb[globdata.indxsampmodi[2+globdata.numbener]] = globdata.modisind
        else:
            globdata.thissampvarb[globdata.indxsampmodi[2:]] = globdata.modispec
            
        
    # death
    if globdata.thisindxprop == globdata.indxpropdeth:
        
        # update the PS index lists
        globdata.thisindxpntsempt[globdata.indxpoplmodi].append(globdata.killindxpnts)
        globdata.thisindxpntsfull[globdata.indxpoplmodi].remove(globdata.killindxpnts)


    # split
    if globdata.thisindxprop == globdata.indxpropsplt:

        # update the PS index lists
        globdata.thisindxpntsfull[globdata.indxpoplmodi].append(globdata.thisindxpntsempt[globdata.indxpoplmodi][0])
        del globdata.thisindxpntsempt[globdata.indxpoplmodi][0]
        
        # update the components
        # first component
        globdata.thissampvarb[globdata.indxsampchd0] = globdata.modilgal[1]
        globdata.thissampvarb[globdata.indxsampchd0+1] = globdata.modibgal[1]
        globdata.thissampvarb[globdata.indxsampchd0+2:globdata.indxsampchd0+2+globdata.numbener] = globdata.modispec[:, 1]
  
        # second component
        globdata.thissampvarb[globdata.indxsampchd1] = globdata.modilgal[2]
        globdata.thissampvarb[globdata.indxsampchd1+1] = globdata.modibgal[2]
        globdata.thissampvarb[globdata.indxsampchd1+2:globdata.indxsampchd1+2+globdata.numbener] = globdata.modispec[:, 2]
        
    # merge
    if globdata.thisindxprop == globdata.indxpropmerg:
        
        # update the PS index lists
        globdata.thisindxpntsfull[globdata.indxpoplmodi].remove(globdata.mergindxchd1)
        globdata.thisindxpntsempt[globdata.indxpoplmodi].append(globdata.mergindxchd1)

        # update the component
        globdata.thissampvarb[globdata.indxsampmodi[0]] = globdata.modilgal[2]
        globdata.thissampvarb[globdata.indxsampmodi[1]] = globdata.modibgal[2]
        globdata.thissampvarb[globdata.indxsampmodi[2:]] = globdata.modispec[:, 2]
        
        
    # component change
    if globdata.thisindxprop >= globdata.indxproplgal:  
        if globdata.thisindxprop == globdata.indxproplgal:
            globdata.thissampvarb[globdata.indxsampmodi] = globdata.modilgal[1]
        elif globdata.thisindxprop == globdata.indxpropbgal:
            globdata.thissampvarb[globdata.indxsampmodi] = globdata.modibgal[1]
        else:
            if globdata.colrprio:
                globdata.thissampvarb[globdata.indxsampmodispec] = globdata.modispec[:, 1]
                if globdata.thisindxprop == globdata.indxpropsind:
                    globdata.thissampvarb[globdata.indxsampmodi] = globdata.modisind
            else:
                globdata.thissampvarb[globdata.indxsampmodi] = globdata.modispec[0, 1]

    # update the unit sample vector
    if globdata.indxsampmodi.size > 0:
        globdata.drmcsamp[globdata.indxsampmodi, 0] = globdata.drmcsamp[globdata.indxsampmodi, 1]


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


def retr_fgl3(globdata):
        
    path = os.environ["PCAT_DATA_PATH"] + '/gll_psc_v16.fit'

    fgl3 = pf.getdata(path)
    
    globdata.fgl3numbpnts = fgl3['glon'].size
    
    globdata.fgl3lgal = fgl3['glon']
    globdata.fgl3lgal = ((globdata.fgl3lgal - 180.) % 360.) - 180.

    globdata.fgl3bgal = fgl3['glat']

    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'] + fgl3['Conf_68_SemiMajor']) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng']) # [rad]
    globdata.fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    globdata.fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    globdata.fgl3sind = fgl3['Spectral_Index']
    
    globdata.fgl3spectype = fgl3['SpectrumType']
    globdata.fgl3scur = fgl3['beta']
    globdata.fgl3scut = fgl3['Cutoff'] * 1e-3
    
    globdata.fgl3timevari = fgl3['Variability_Index']
    globdata.indxfgl3timevari = where(globdata.fgl3timevari > 72.44)[0]
    
    fgl3spectemp = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], fgl3['Flux10000_100000']))
    fgl3spectemp = fgl3spectemp[globdata.indxenerincl, :] / globdata.diffener[:, None]
    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], fgl3['Unc_Flux10000_100000']))
    fgl3specstdvtemp = fgl3specstdvtemp[globdata.indxenerincl, :, :] / globdata.diffener[:, None, None]
    
    globdata.fgl3spec = zeros((3, globdata.numbener, globdata.fgl3numbpnts))
    globdata.fgl3spec[0, :, :] = fgl3spectemp
    globdata.fgl3spec[1, :, :] = fgl3spectemp - fgl3specstdvtemp[:, :, 0]
    globdata.fgl3spec[2, :, :] = fgl3spectemp + fgl3specstdvtemp[:, :, 1]
    
    indxpixlfgl3 = retr_indxpixl(globdata, globdata.fgl3bgal, globdata.fgl3lgal)
    globdata.fgl3cnts = globdata.fgl3spec[0, :, :, None] * globdata.expo[:, indxpixlfgl3, :] * globdata.diffener[:, None, None]
    globdata.fgl3gang = rad2deg(arccos(cos(deg2rad(globdata.fgl3lgal)) * cos(deg2rad(globdata.fgl3bgal))))
        
    # adjust 3FGL positions according to the ROI center
    if globdata.regitype == 'ngal':
        rttr = hp.rotator.Rotator(rot=[0., 90., 0.], deg=True)
        globdata.fgl3bgal, globdata.fgl3lgal = rad2deg(rttr(deg2rad(90. - globdata.fgl3bgal), deg2rad(globdata.fgl3lgal)))
        globdata.fgl3bgal = 90. - globdata.fgl3bgal

    # select the 3FGL point sources in the ROI
    globdata.indxfgl3rofi = arange(globdata.fgl3lgal.size, dtype=int)
    for i in globdata.indxener:
        globdata.indxfgl3rofi = intersect1d(where((globdata.fgl3spec[0, i, :] > globdata.minmspec[i]) & \
            (globdata.fgl3spec[0, i, :] < globdata.maxmspec[i]))[0], globdata.indxfgl3rofi)
    globdata.indxfgl3rofi = intersect1d(where((abs(globdata.fgl3lgal) < globdata.maxmgangmarg) & \
            (abs(globdata.fgl3bgal) < globdata.maxmgangmarg))[0], globdata.indxfgl3rofi)
    globdata.fgl3numbpntsrofi = globdata.indxfgl3rofi.size
    globdata.indxfgl3timevarirofi = where(globdata.fgl3timevari[globdata.indxfgl3rofi] > 72.44)[0]

    # sanity check
    for i in globdata.indxener:
        if (globdata.fgl3spec[0, i, globdata.indxfgl3rofi] > globdata.maxmspec[i]).any():
            print 'maxmspec %d is bad!' % i


def retr_rtag(globdata, indxprocwork):
    
    if indxprocwork == None:
        rtag = 'AA_%d_%d_%d_%d_%s_%s_%s' % (globdata.numbproc, globdata.numbswep, globdata.numbburn, globdata.factthin, \
            globdata.datatype, globdata.regitype, globdata.psfntype)
    else:
        rtag = '%02d_%d_%d_%d_%d_%s_%s_%s' % (indxprocwork, globdata.numbproc, globdata.numbswep, globdata.numbburn, globdata.factthin, \
            globdata.datatype, globdata.regitype, globdata.psfntype)
        
    return rtag


def retr_gaus(globdata, indxsamp, stdv):
    
    if globdata.fracrand > 0.:
        if rand() < globdata.fracrand:
            globdata.drmcsamp[indxsamp, 1] = rand()
        else:
            globdata.drmcsamp[indxsamp, 1] = globdata.drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        globdata.drmcsamp[indxsamp, 1] = globdata.drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_dist(globdata, lgal0, bgal0, lgal1, bgal1):
    
    if globdata.pixltype == 'heal':
        dir1 = array([lgal0, bgal0])
        dir2 = array([lgal1, bgal1])
        dist = angdist(dir1, dir2, lonlat=True)
    else:
        dist = deg2rad(sqrt((lgal0 - lgal1)**2 + (bgal0 - bgal1)**2))

    return dist
       
    
def retr_strgprop(globdata):
    
    globdata.strgprop = ['FDF Norm', 'FDF Slope', 'PSF', 'Back. Norm', 'Birth', 'Death', 'Split', 'Merge', \
  'Longitude Update', 'Latitude Update', 'Flux Update', 'Spectral Index Update']


def retr_singgaus(scaldevi, sigc):
    
    psfn = 1. / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_singking(scaldevi, sigc, gamc):
    
    psfn = 1. / 2. / pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc)

    return psfn


def retr_doubgaus(scaldevi, frac, sigc, sigt):
    
    psfn = frac / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2) +                 (1. - frac) / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2)

    return psfn


def retr_gausking(scaldevi, frac, sigc, sigt, gamt):

    psfn = frac / 2. / pi / sigc**2 * exp(-0.5 * scaldevi**2 / sigc**2) +                 (1. - frac) / 2. / pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_doubking(scaldevi, frac, sigc, gamc, sigt, gamt):

    psfn = frac / 2. / pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc) +                 (1. - frac) / 2. / pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def chsq_fdfnslop(globdata, para, i):

    fdfnslop = para[0]
    fdfnnormnorm = para[1]
    
    fluxhistmodl = fdfnnormnorm * globdata.diffspec[i, :] * pdfn_spec(globdata.meanspec[i, :], fdfnslop, globdata.minmspec[i], globdata.maxmspec[i])

    chsq = sum(((fluxhistmodl.flatten()[jspecfgl3] - globdata.fgl3spechist[i, jspecfgl3]) / globdata.fgl3spechist[i, jspecfgl3])**2)
    
    return chsq


def retr_strgfluxunit(globdata):
    
    if globdata.exprtype == 'ferm':
        strgfluxunit = r'[1/cm$^2$/s/GeV]'
    if globdata.exprtype == 'sdss':
        strgfluxunit = '[nMgy]'
        
    return strgfluxunit
     
    
def retr_enerstrg(globdata):
    
    if globdata.exprtype == 'ferm':
        binsenerstrg = []
        enerstrg = []
        for i in globdata.indxener:
            binsenerstrg.append('%.3g GeV - %.3g GeV' % (globdata.binsener[i], globdata.binsener[i+1]))
            enerstrg.append('%.3g GeV' % globdata.meanener[i])
            
    if globdata.exprtype == 'sdss':
        globdata.binsenerstrg = ['i-band', 'r-band', 'g-band']
        enerstrg = ['i-band', 'r-band', 'g-band']
        
    return enerstrg, binsenerstrg


def retr_prop(globdata):
  
    globdata.thisindxsamplgal, globdata.thisindxsampbgal,         globdata.thisindxsampspec, globdata.thisindxsampsind,         globdata.thisindxsampcomp = retr_indx(globdata, globdata.thisindxpntsfull)
    
    if globdata.verbtype > 2:
        print 'retr_prop(): '

        print 'drmcsamp'
        print globdata.drmcsamp
        
        print 'thissampvarb: '
        for k in range(globdata.thissampvarb.size):
            if k == globdata.indxsampcompinit:
                print
            if k > globdata.indxsampcompinit and (k - globdata.indxsampcompinit) % globdata.numbcomp == 0:
                print
            print globdata.thissampvarb[k]
        print
            
        print 'thisindxpntsfull: ', globdata.thisindxpntsfull
        print 'thisindxpntsempt: ', globdata.thisindxpntsempt  
        print 'thisindxsamplgal: ', globdata.thisindxsamplgal
        print 'thisindxsampbgal: ', globdata.thisindxsampbgal
        print 'thisindxsampspec: '
        print globdata.thisindxsampspec
        if globdata.colrprio:
            print 'thisindxsampsind: ', globdata.thisindxsampsind
        print 'thisindxsampcomp: ', globdata.thisindxsampcomp
        print
        
    # hyper-parameter changes
    # mean number of point sources
    if globdata.thisindxprop == globdata.indxpropfdfnnorm:
        globdata.indxsampmodi = globdata.indxsampfdfnnorm[globdata.indxpoplmodi]
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvfdfnnorm)
        globdata.nextsampvarb[globdata.indxsampfdfnnorm] = icdf_logt(globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.minmfdfnnorm, globdata.factfdfnnorm)
        if globdata.colrprio:
            globdata.indxenermodi = globdata.indxenerfdfn
        else:
            globdata.indxenermodi = globdata.indxener
        
    # flux distribution function slope
    if globdata.thisindxprop == globdata.indxpropfdfnslop:
        if globdata.colrprio:
            globdata.indxenermodi = globdata.indxenerfdfn
        else:
            globdata.indxenermodi = choice(globdata.indxener)
        globdata.indxsampmodi = globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvfdfnslop)
        
        globdata.nextsampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]] = \
            icdf_atan(globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.minmfdfnslop, globdata.factfdfnslop)
        if globdata.colrprio:
            globdata.indxsampmodi = concatenate((globdata.indxsampmodi, globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxener, :].flatten()))
        else:
            globdata.indxsampmodi = concatenate((array([globdata.indxsampmodi]), globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, :]))
            
            
        if globdata.verbtype > 2:
            print 'indxpoplmodi'
            print globdata.indxpoplmodi
            print 'indxenermodi'
            print globdata.indxenermodi
            print 'nextsampvarb[globdata.indxsampfdfnslop]'
            print globdata.nextsampvarb[globdata.indxsampfdfnslop]
            print 'indxsampmodi'
            print globdata.indxsampmodi
        
        
            
    # PSF parameter change 
    if globdata.thisindxprop == globdata.indxproppsfipara:
        
        # index of the PSF parameter to change
        globdata.indxpsfiparamodi = choice(globdata.indxpsfipara)

        # the energy bin of the PS flux map to be modified
        globdata.indxenermodi = array([(globdata.indxpsfiparamodi % globdata.numbpsfiparaevtt) // globdata.numbformpara])
        
        # sample index to be modified
        globdata.indxsampmodi = globdata.indxsamppsfipara[globdata.indxpsfiparamodi]
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvpsfipara)
        globdata.nextsampvarb[globdata.indxsamppsfipara] = copy(globdata.thissampvarb[globdata.indxsamppsfipara])
        globdata.nextsampvarb[globdata.indxsamppsfipara[globdata.indxpsfiparamodi]] = \
            icdf_psfipara(globdata, globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.indxpsfiparamodi)
        globdata.modilgal = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]]
        globdata.modibgal = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]]
        globdata.modispec = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi]]
        
        if globdata.verbtype > 1:
           
            print 'indxpsfiparamodi: ', globdata.indxpsfiparamodi
            print 'indxenermodi: ', globdata.indxenermodi
            print 'indxsampmodi: ', globdata.indxsampmodi
            print 'globdata.drmcsamp[globdata.indxsampmodi, -1]'
            print globdata.drmcsamp[globdata.indxsampmodi, -1]
            print 'icdf_psfipara(globdata, globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.indxpsfiparamodi)'
            print icdf_psfipara(globdata, globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.indxpsfiparamodi)
            print 'thispsfipara: ', globdata.thissampvarb[globdata.indxsamppsfipara]
            print 'thissampvarb[indxsampmodi]: ', globdata.thissampvarb[globdata.indxsampmodi]
            print 'nextpsfipara: ', globdata.nextsampvarb[globdata.indxsamppsfipara]
            print 'nextpsfipara[indxsampmodi]: ', globdata.nextsampvarb[globdata.indxsampmodi]
            print 

        
    # background changes
    
    # diffuse model
    if globdata.thisindxprop == globdata.indxpropnormback:

        # determine the sample index to be changed
        globdata.indxenermodi = choice(globdata.indxener)
        globdata.indxbackmodi = choice(globdata.indxback)
        globdata.indxsampmodi = globdata.indxsampnormback[globdata.indxbackmodi, globdata.indxenermodi]
        
        # propose
        retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvback)

        # transform back from the unit space
        globdata.nextsampvarb[globdata.indxsampmodi] = icdf_logt(globdata.drmcsamp[globdata.indxsampmodi, -1], \
            globdata.minmnormback[globdata.indxbackmodi], globdata.factnormback[globdata.indxbackmodi])

        if globdata.verbtype > 1:
            print 'indxsampmodi: ', globdata.indxsampmodi
            print 'nextsampvarb[globdata.indxsampmodi]: ', globdata.nextsampvarb[globdata.indxsampmodi]

    
    # birth
    if globdata.thisindxprop == globdata.indxpropbrth:

        # change the number of PS
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] + 1
    
        # initial sample index to add the new PS
        indxsampbrth = int(globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + globdata.thisindxpntsempt[globdata.indxpoplmodi][0] * globdata.numbcomp)
        
        # sample auxiliary variables
        if globdata.colrprio:
            numbauxipara = globdata.numbcompcolr
        else:
            numbauxipara = globdata.numbcomp
        globdata.auxipara = rand(numbauxipara)

        if globdata.colrprio:
            globdata.drmcsamp[indxsampbrth:indxsampbrth+2, -1] = globdata.auxipara[0:2]
            globdata.drmcsamp[indxsampbrth+2+globdata.indxenerfdfn, -1] = globdata.auxipara[-2]
            globdata.drmcsamp[indxsampbrth+globdata.numbcomp-1, -1] = globdata.auxipara[-1]
        else:   
            globdata.drmcsamp[indxsampbrth:indxsampbrth+globdata.numbcomp, -1] = globdata.auxipara

        # sample indices to be modified
        globdata.indxsampmodi = arange(indxsampbrth, indxsampbrth + globdata.numbcomp, dtype=int)

        # modification catalog
        globdata.modilgal = empty(1)
        globdata.modibgal = empty(1)
        if globdata.colrprio:
            globdata.modiflux = empty(1)
            globdata.modisind = empty(1)
        globdata.modispec = zeros((globdata.numbener, 1))
        
        globdata.modilgal[0] = icdf_self(globdata.drmcsamp[indxsampbrth, -1], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
        globdata.modibgal[0] = icdf_self(globdata.drmcsamp[indxsampbrth+1, -1], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)

        if globdata.colrprio:
            globdata.modiflux[0] = icdf_spec(globdata, globdata.drmcsamp[indxsampbrth+2+globdata.indxenerfdfn, -1], \
                    globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenerfdfn]], \
                    globdata.minmspec[globdata.indxenerfdfn], globdata.maxmspec[globdata.indxenerfdfn])
            globdata.modisind[0] = icdf_atan(globdata.drmcsamp[indxsampbrth+globdata.numbcomp-1, -1], globdata.minmsind, globdata.factsind)
            globdata.modispec[:, 0] = retr_spec(globdata, globdata.modiflux[0], globdata.modisind[0]).flatten()
        else:
            for i in globdata.indxener:
                globdata.modispec[i, 0] = icdf_spec(globdata, globdata.drmcsamp[indxsampbrth+2+i, -1], \
                    globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])
    
        if globdata.verbtype > 1:
            print 'auxipara: ', globdata.auxipara
            print 'modilgal: ', globdata.modilgal
            print 'modibgal: ', globdata.modibgal
            print 'modispec: '
            print globdata.modispec
            print
            
    # kill
    if globdata.thisindxprop == globdata.indxpropdeth:
        
        # change the number of PS
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] - 1

        # occupied PS index to be killed
        killindxindxpnts = choice(arange(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        globdata.killindxpnts = globdata.thisindxpntsfull[globdata.indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        globdata.indxsampmodi = array([])
            
        # modification catalog
        globdata.modilgal = empty(1)
        globdata.modibgal = empty(1)
        globdata.modispec = zeros((globdata.numbener, 1))
        globdata.modilgal[0] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][killindxindxpnts]]
        globdata.modibgal[0] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][killindxindxpnts]]
        globdata.modispec[:, 0] = -globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, killindxindxpnts]]

        if globdata.verbtype > 1:
            print 'killindxpnts: ', globdata.killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', globdata.modilgal
            print 'modibgal: ', globdata.modibgal
            print 'modispec: '
            print globdata.modispec
            print
            
  
    # split
    if globdata.thisindxprop == globdata.indxpropsplt:
        
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] + 1
        
        # determine which point source to split
        spltindxindxpnts = choice(arange(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]], dtype=int))
        spltindxpnts = globdata.thisindxpntsfull[globdata.indxpoplmodi][spltindxindxpnts]
    
        # update the sample vector
        globdata.indxsampchd0 = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + globdata.thisindxpntsfull[globdata.indxpoplmodi][spltindxindxpnts] * globdata.numbcomp
        indxfinl0 = globdata.indxsampchd0 + globdata.numbcomp
        globdata.indxsampchd1 = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + globdata.thisindxpntsempt[globdata.indxpoplmodi][0] * globdata.numbcomp
        indxfinl1 = globdata.indxsampchd1 + globdata.numbcomp
        
        # determine the modified sample vector indices
        globdata.indxsampmodi = concatenate((arange(globdata.indxsampchd0, indxfinl0, dtype=int), arange(globdata.indxsampchd1, indxfinl1, dtype=int)))
        
        thislgal = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][spltindxindxpnts]]
        thisbgal = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][spltindxindxpnts]]
        thisspec = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, spltindxindxpnts]]
        
        if globdata.verbtype > 1:
            print 'spltindxindxpnts: ', spltindxindxpnts
            print 'spltindxpnts: ', spltindxpnts
            print 'indxsampchd0: ', globdata.indxsampchd0
            print 'indxfinl0: ', indxfinl0
            print 'indxsampchd1: ', globdata.indxsampchd1
            print 'indxfinl1: ', indxfinl1
            if pixltype == 'heal':
                print 'thislgal: ', thislgal
                print 'thisbgal: ', thisbgal
            else:
                print 'thislgal: ', 3600. * thislgal
                print 'thisbgal: ', 3600. * thisbgal
            print 'thisspec: ', thisspec
            
            
        # determine the new components
        globdata.auxipara = empty(globdata.numbcomp)
        globdata.auxipara[0:2] = rand(2) * spmrlbhl
        globdata.auxipara[2:] = (exp(rand(globdata.numbener)) - 1.) / (exp(1.) - 1.) * (globdata.maxmspec - globdata.minmspec) + globdata.minmspec
        
        if globdata.verbtype > 1:
            if pixltype == 'heal':
                print 'auxipara[0]: ', globdata.auxipara[0]
                print 'auxipara[1]: ', globdata.auxipara[1]
            else:
                print 'auxipara[0]: ', 3600. * globdata.auxipara[0]
                print 'auxipara[1]: ', 3600. * globdata.auxipara[1]
            print 'auxipara[2:]: ', globdata.auxipara[2:]
            print
            
        nextlgal0 = thislgal + globdata.auxipara[0]
        nextlgal1 = thislgal - globdata.auxipara[0]
        nextbgal0 = thisbgal + globdata.auxipara[1]
        nextbgal1 = thisbgal - globdata.auxipara[1]
        nextspec0 = (thisspec + globdata.auxipara[2:]) / 2.
        nextspec1 = (thisspec - globdata.auxipara[2:]) / 2.
        
        if globdata.verbtype > 1:
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

        if abs(nextlgal0) > globdata.maxmgangmarg or abs(nextlgal1) > globdata.maxmgangmarg or \
           abs(nextbgal0) > globdata.maxmgangmarg or abs(nextbgal1) > globdata.maxmgangmarg or \
           where((nextspec0 > globdata.maxmspec) | (nextspec0 < globdata.minmspec))[0].size > 0 or \
           where((nextspec1 > globdata.maxmspec) | (nextspec1 < globdata.minmspec))[0].size > 0:
               globdata.reje = True
                
        if not globdata.reje:

            lgal = concatenate((array([nextlgal0, nextlgal1]), setdiff1d(globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([nextbgal0, nextbgal1]), setdiff1d(globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]], thisbgal)))
            pairlist = retr_pairlist(lgal, bgal)

            ## first new component
            globdata.drmcsamp[globdata.indxsampchd0, -1] = cdfn_self(nextlgal0, -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            globdata.drmcsamp[globdata.indxsampchd0+1, -1] = cdfn_self(nextbgal0, -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            for i in globdata.indxener:
                globdata.drmcsamp[globdata.indxsampchd0+2+i, -1] = cdfn_spec(nextspec0[i], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])

            ## second new component
            globdata.drmcsamp[globdata.indxsampchd1, -1] = cdfn_self(nextlgal1, -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            globdata.drmcsamp[globdata.indxsampchd1+1, -1] = cdfn_self(nextbgal1, -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            for i in globdata.indxener:
                globdata.drmcsamp[globdata.indxsampchd1+2+i, -1] = cdfn_spec(nextspec1[i], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])


            # construct the modification catalog
            globdata.modilgal = empty(3)
            globdata.modibgal = empty(3)
            globdata.modispec = empty((globdata.numbener, 3))

            ## component to be removed
            globdata.modilgal[0] = thislgal
            globdata.modibgal[0] = thisbgal
            globdata.modispec[:, 0] = -thisspec.flatten()

            ## first component to be added
            globdata.modilgal[1] = nextlgal0
            globdata.modibgal[1] = nextbgal0
            globdata.modispec[:, 1] = nextspec0.flatten()

            # second component to be added
            globdata.modilgal[2] = nextlgal1
            globdata.modibgal[2] = nextbgal1
            globdata.modispec[:, 2] = nextspec1.flatten()

    if globdata.thisindxprop == globdata.indxpropmerg:
        
        globdata.nextsampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] = globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]] - 1

        # determine the first PS to merge
        #dir2 = array([globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]], globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]]])
            
        lgal = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi]]
        bgal = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi]]
        pairlist = retr_pairlist(lgal, bgal)
        
        if globdata.verbtype > 1:
            print 'lgal'
            print lgal
            print 'bgal'
            print bgal
            print 'pairlist'
            print pairlist
            
            
        if len(pairlist) == 0:
            globdata.reje = True
        else:
            globdata.reje = False
            jpair = choice(arange(len(pairlist)))
            mergindxindxpnts0 = pairlist[jpair][0]
            mergindxindxpnts1 = pairlist[jpair][1]
  
        if not globdata.reje:

            # fisrt PS index to be merged
            mergindxchd0 = globdata.thisindxpntsfull[globdata.indxpoplmodi][mergindxindxpnts0]
            mergindxsampinit0 = globdata.indxsampcompinit + mergindxchd0 * globdata.numbcomp

            # second PS index to be merged
            globdata.mergindxchd1 = globdata.thisindxpntsfull[globdata.indxpoplmodi][mergindxindxpnts1]
            mergindxsampinit1 = globdata.indxsampcompinit + globdata.mergindxchd1 * globdata.numbcomp

            # determine the modified sample vector indices
            globdata.indxsampchd0 = globdata.indxsampcompinit + globdata.numbcomp * mergindxchd0
            indxfinl0 = globdata.indxsampchd0 + globdata.numbcomp
            globdata.indxsampchd1 = globdata.indxsampcompinit + globdata.numbcomp * globdata.mergindxchd1
            indxfinl1 = globdata.indxsampchd1 + globdata.numbcomp

            globdata.indxsampmodi = arange(globdata.indxsampchd0, indxfinl0)

            # indices of the PS to be merges
            mergindxpnts = sort(array([mergindxchd0, globdata.mergindxchd1], dtype=int))

            thislgal0 = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][mergindxindxpnts0]]
            thisbgal0 = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][mergindxindxpnts0]]
            thisspec0 = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, mergindxindxpnts0]]

            thislgal1 = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][mergindxindxpnts1]]
            thisbgal1 = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][mergindxindxpnts1]]
            thisspec1 = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][:, mergindxindxpnts1]]

            # auxiliary component
            globdata.auxipara = zeros(globdata.numbcomp)
            globdata.auxipara[0] = (thislgal0 - thislgal1) / 2.
            globdata.auxipara[1] = (thisbgal0 - thisbgal1) / 2.
            globdata.auxipara[2:] = thisspec0 - thisspec1

            # merged PS
            nextlgal = (thislgal0 + thislgal1) / 2.
            nextbgal = (thisbgal0 + thisbgal1) / 2.
            nextspec = thisspec0 + thisspec1
            
            globdata.drmcsamp[globdata.indxsampchd0, -1] = cdfn_self(nextlgal, -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            globdata.drmcsamp[globdata.indxsampchd0+1, -1] = cdfn_self(nextbgal, -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            for i in globdata.indxener:
                globdata.drmcsamp[globdata.indxsampchd0+2+i, -1] = cdfn_spec(nextspec[i], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, i]], globdata.minmspec[i], globdata.maxmspec[i])

            # construct the modification catalog
            globdata.modilgal = empty(3)
            globdata.modibgal = empty(3)
            globdata.modispec = empty((globdata.numbener, 3))

            ## first component to be merged
            globdata.modilgal[0] = thislgal0
            globdata.modibgal[0] = thisbgal0
            globdata.modispec[:, 0] = -thisspec0.flatten()

            ## first component to be merged
            globdata.modilgal[1] = thislgal1
            globdata.modibgal[1] = thisbgal1
            globdata.modispec[:, 1] = -thisspec1.flatten()

            ## component to be added
            globdata.modilgal[2] = nextlgal
            globdata.modibgal[2] = nextbgal
            globdata.modispec[:, 2] = nextspec.flatten()

            if globdata.verbtype > 1:
                print 'mergindxchd0: ', mergindxchd0
                print 'mergindxindxpnts0: ', mergindxindxpnts0
                print 'mergindxchd1: ', globdata.mergindxchd1
                print 'mergindxindxpnts1: ', mergindxindxpnts1
                print 'indxsampchd0: ', globdata.indxsampchd0
                print 'indxfinl0: ', indxfinl0
                print 'indxsampchd1: ', globdata.indxsampchd1
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
                    print 'auxipara[0]: ', globdata.auxipara[0]
                    print 'auxipara[1]: ', globdata.auxipara[1]
                else:
                    print 'nextlgal: ', 3600. * nextlgal
                    print 'nextbgal: ', 3600. * nextbgal
                    print 'auxipara[0]: ', 3600. * globdata.auxipara[0]
                    print 'auxipara[1]: ', 3600. * globdata.auxipara[1]
                print 'nextspec: ', nextspec
                print 'auxipara[2:]: ', globdata.auxipara[2:]
                print

    # component change
    if globdata.thisindxprop >= globdata.indxproplgal:     
        
        if globdata.thisindxprop == globdata.indxproplgal or globdata.thisindxprop == globdata.indxpropbgal:
            if globdata.thisindxprop == globdata.indxproplgal:
                globdata.indxcompmodi = 0
            else:
                globdata.indxcompmodi = 1
            globdata.indxenermodi = globdata.indxener
        else:
            if globdata.colrprio:
                globdata.indxenermodi = globdata.indxener
                if globdata.thisindxprop == globdata.indxpropspec:
                    globdata.indxcompmodi = 2 + globdata.indxenerfdfn
                elif globdata.thisindxprop == globdata.indxpropsind:
                    globdata.indxcompmodi = array([2 + globdata.numbener])
            else:
                globdata.indxenermodi = choice(globdata.indxener)
                globdata.indxcompmodi = globdata.indxenermodi + 2
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = globdata.thisindxpntsfull[globdata.indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        globdata.indxsampmodiinit = globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.indxpoplmodi + modiindxpnts * globdata.numbcomp
        
        # sample index to be modified
        globdata.indxsampmodi = globdata.indxsampmodiinit + globdata.indxcompmodi
        if globdata.colrprio:
            globdata.indxsampmodispec = globdata.indxsampmodiinit + 2 + globdata.indxener
        
        # propose
        if globdata.thisindxprop == globdata.indxpropspec:
            retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvspec)
        else:
            retr_gaus(globdata, globdata.indxsampmodi, globdata.stdvlbhl) 

        # modification catalog
        globdata.modilgal = empty(2)
        globdata.modibgal = empty(2)
        globdata.modispec = empty((globdata.indxenermodi.size, 2))
  
        if globdata.colrprio:
            thisflux = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenerfdfn, modiindxindxpnts]]
            thissind = globdata.thissampvarb[globdata.thisindxsampsind[globdata.indxpoplmodi][modiindxindxpnts]]
            thisspec = retr_spec(globdata, thisflux, thissind)
        else:
            thisspec = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenermodi, modiindxindxpnts]]
            
        globdata.modispec[:, 0] = -thisspec.flatten()
        if globdata.thisindxprop == globdata.indxproplgal or globdata.thisindxprop == globdata.indxpropbgal:
            if globdata.indxcompmodi == 0:
                globdata.modilgal[0] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][modiindxindxpnts]]
                globdata.modilgal[1] = icdf_self(globdata.drmcsamp[globdata.indxsampmodi, -1], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
                globdata.modibgal[:] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][modiindxindxpnts]]
            else:
                globdata.modilgal[:] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][modiindxindxpnts]]
                globdata.modibgal[0] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][modiindxindxpnts]]
                globdata.modibgal[1] = icdf_self(globdata.drmcsamp[globdata.indxsampmodi, -1], -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            globdata.modispec[:, 1] = thisspec.flatten()
        else:
            globdata.modilgal[:] = globdata.thissampvarb[globdata.thisindxsamplgal[globdata.indxpoplmodi][modiindxindxpnts]]
            globdata.modibgal[:] = globdata.thissampvarb[globdata.thisindxsampbgal[globdata.indxpoplmodi][modiindxindxpnts]]
            if globdata.colrprio:
                if globdata.thisindxprop == globdata.indxpropspec:
                    globdata.modiflux = icdf_spec(globdata, globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenerfdfn]], 
               globdata.minmspec[globdata.indxenerfdfn], globdata.maxmspec[globdata.indxenerfdfn])
                    globdata.modisind = globdata.thissampvarb[globdata.thisindxsampsind[globdata.indxpoplmodi][modiindxindxpnts]]
                else:
                    globdata.modiflux = globdata.thissampvarb[globdata.thisindxsampspec[globdata.indxpoplmodi][globdata.indxenerfdfn, modiindxindxpnts]]
                    globdata.modisind = icdf_atan(globdata.drmcsamp[globdata.indxsampmodi, -1], globdata.minmsind, globdata.factsind)

                globdata.modispec[:, 1] = retr_spec(globdata, globdata.modiflux, globdata.modisind).flatten()
            else:
                specunit = globdata.drmcsamp[globdata.indxsampmodi, -1]
                fdfnslop = globdata.thissampvarb[globdata.indxsampfdfnslop[globdata.indxpoplmodi, globdata.indxenermodi]]
                globdata.modispec[:, 1] = icdf_spec(globdata, specunit, fdfnslop, globdata.minmspec[globdata.indxenermodi], globdata.maxmspec[globdata.indxenermodi])

        # log
        if globdata.verbtype > 1:
            print 'modilgal: ', globdata.modilgal
            print 'modibgal: ', globdata.modibgal
            print 'modispec: '
            print globdata.modispec
            print 'indxcompmodi: ', globdata.indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts


    # energy bin in which to evaluate the log-likelihood
    if globdata.indxpropbrth <= globdata.thisindxprop <= globdata.indxpropmerg:
        globdata.indxenermodi = arange(globdata.numbener)

    if type(globdata.indxenermodi) == int64:
        globdata.indxenermodi = array([globdata.indxenermodi])

    if globdata.verbtype > 1:
        print 'indxsampmodi: ', globdata.indxsampmodi
        print 'indxenermodi: ', globdata.indxenermodi

    # auxiliary variable density fraction and jacobian
    if (globdata.thisindxprop == globdata.indxpropsplt or globdata.thisindxprop == globdata.indxpropmerg) and not globdata.reje:

        spltglobdata.combfact = log(globdata.thissampvarb[globdata.indxsampnumbpnts[globdata.indxpoplmodi]]**2 / len(pairlist))
        
        if globdata.thisindxprop == globdata.indxpropsplt:
            thisglobdata.combfact = spltglobdata.combfact 
            thisglobdata.jcbnfact = spltglobdata.jcbnfact
        else:
            thisglobdata.combfact = -spltglobdata.combfact 
            thisglobdata.jcbnfact = -spltglobdata.jcbnfact


        globdata.laccfrac = thisglobdata.jcbnfact + thisglobdata.combfact

        globdata.listglobdata.numbpair[globdata.thisindxswep] = len(pairlist)
        globdata.listglobdata.jcbnfact[globdata.thisindxswep] = thisglobdata.jcbnfact
        globdata.listglobdata.combfact[globdata.thisindxswep] = thisglobdata.combfact
        globdata.listglobdata.auxipara[globdata.thisindxswep, :] = globdata.auxipara
        globdata.listglobdata.laccfrac[globdata.thisindxswep] = globdata.laccfrac

    else:
        globdata.laccfrac = 0.  
        
        
def retr_psfn(globdata, psfipara, indxenertemp, thisangl, psfntype):

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

    indxpsfiparatemp = numbformpara * (indxenertemp[:, None] + globdata.numbener * globdata.indxevtt[None, :])
    
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


def retr_psfimodl(globdata):
    
    # PSF parameters
    if globdata.psfntype == 'singgaus':
        globdata.numbformpara = 1
    elif globdata.psfntype == 'singking':
        globdata.numbformpara = 2 
    elif globdata.psfntype == 'doubgaus':
        globdata.numbformpara = 3
    elif globdata.psfntype == 'gausking':
        globdata.numbformpara = 4
    elif globdata.psfntype == 'doubking':
        globdata.numbformpara = 5
       
    globdata.indxformpara = arange(globdata.numbformpara) 
    globdata.numbpsfiparaevtt = globdata.numbener * globdata.numbformpara
    globdata.numbpsfipara = globdata.numbpsfiparaevtt * globdata.numbevtt
    globdata.indxpsfipara = arange(globdata.numbpsfipara)
    globdata.indxmodlpsfipara = arange(globdata.numbpsfipara)   

    minmformpara = zeros(globdata.numbformpara)
    maxmformpara = zeros(globdata.numbformpara)
    factformpara = zeros(globdata.numbformpara)
    scalformpara = zeros(globdata.numbformpara, dtype=object)
    if globdata.exprtype == 'ferm':
        minmanglpsfn = deg2rad(0.01)
        maxmanglpsfn = deg2rad(3.)
        minmgamm = 2.
        maxmgamm = 20.
        minmpsfnfrac = 0.
        maxmpsfnfrac = 1.
    if globdata.exprtype == 'sdss':
        minmanglpsfn = deg2rad(0.01 / 3600.)
        maxmanglpsfn = deg2rad(2. / 3600.)
    if globdata.psfntype == 'singgaus':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        scalformpara[0] = 'logt'
        strgformpara = [r'$\sigma']
    elif globdata.psfntype == 'singking':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        minmformpara[1] = minmgamm
        maxmformpara[1] = maxmgamm
        scalformpara[0] = 'logt'
        scalformpara[1] = 'atan'
        strgformpara = [r'$\sigma', r'$\gamma']
    elif globdata.psfntype == 'doubgaus':
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
    elif globdata.psfntype == 'gausking':
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
    elif globdata.psfntype == 'doubking':
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

    for k in globdata.indxformpara:
        if scalformpara[k] == 'self':
            factformpara[k] = maxmformpara[k] - minmformpara[k]
        if scalformpara[k] == 'logt':
            factformpara[k] = log(maxmformpara[k] / minmformpara[k])
        if scalformpara[k] == 'atan':
            factformpara[k] = arctan(maxmformpara[k]) - arctan(minmformpara[k])
            
    globdata.minmpsfipara = tile(tile(minmformpara, globdata.numbener), globdata.numbevtt)
    globdata.maxmpsfipara = tile(tile(maxmformpara, globdata.numbener), globdata.numbevtt)
    globdata.scalpsfipara = tile(tile(scalformpara, globdata.numbener), globdata.numbevtt)
    globdata.factpsfipara = tile(tile(factformpara, globdata.numbener), globdata.numbevtt)
    
    # PSF parameter strings
    globdata.strgpsfipara = [strgformpara[k] + '^{%d%d}$' % (globdata.indxenerincl[i], globdata.indxevttincl[m]) \
        for m in globdata.indxevtt for i in globdata.indxener for k in globdata.indxformpara]
    globdata.indxpsfiparainit = (globdata.indxevtt[:, None] * globdata.numbpsfiparaevtt + globdata.indxener[None, :] * globdata.numbformpara).flatten()
    for k in arange(globdata.indxpsfiparainit.size):
        if globdata.psfntype == 'singgaus' or globdata.psfntype == 'singking':
            globdata.strgpsfipara[globdata.indxpsfiparainit[k]] += ' ' + globdata.strganglunit
        elif globdata.psfntype == 'doubgaus' or globdata.psfntype == 'gausking':
            globdata.strgpsfipara[globdata.indxpsfiparainit[k]+1] += ' ' + globdata.strganglunit
            globdata.strgpsfipara[globdata.indxpsfiparainit[k]+2] += ' ' + globdata.strganglunit
        elif globdata.psfntype == 'doubking':
            globdata.strgpsfipara[globdata.indxpsfiparainit[k]+1] += ' ' + globdata.strganglunit
            globdata.strgpsfipara[globdata.indxpsfiparainit[k]+3] += ' ' + globdata.strganglunit
            
    globdata.indxpsfipara = arange(globdata.numbpsfipara)



def retr_propmodl(globdata):

    globdata.strgprop = ['fdfnnorm', 'fdfnslop', 'psfipara', 'normback', 'brth', 'deth', 'splt', 'merg', 'lgal', 'bgal', 'spec', 'sind']

    globdata.numbprop = len(globdata.strgprop)
    globdata.indxprop = arange(globdata.numbprop)
    globdata.indxpropfdfnnorm = 0
    globdata.indxpropfdfnslop = 1
    globdata.indxproppsfipara = 2
    globdata.indxpropnormback = 3
    globdata.indxpropbrth = 4
    globdata.indxpropdeth = 5
    globdata.indxpropsplt = 6
    globdata.indxpropmerg = 7
    globdata.indxproplgal = 8
    globdata.indxpropbgal = 9
    globdata.indxpropspec = 10
    globdata.indxpropsind = 11

    if globdata.probprop == None:

        probfdfnnorm = array([1.])
        if globdata.colrprio:
            probfdfnslop = array([1.])
        else:
            #probfdfnslop = array([1.] * globdata.numbener)
            probfdfnslop = array([1.])
            
        #probpsfipara = array([1.] * globdata.numbpsfipara)
        #probnormback = array([1.] * globdata.numbback * globdata.numbener)
        
        if globdata.proppsfn:
            probpsfipara = array([1.])
        else:
            probpsfipara = array([0.])
        probnormback = array([1.])
        
        probbrth = array([0.1 * sum(globdata.maxmnumbpnts)])
        probdeth = array([0.1 * sum(globdata.maxmnumbpnts)])
        probsplt = array([0. * sum(globdata.maxmnumbpnts)])
        probmerg = array([0. * sum(globdata.maxmnumbpnts)])
        
        problgal = array([sum(globdata.maxmnumbpnts) / 2.])
        probbgal = array([sum(globdata.maxmnumbpnts) / 2.])
        if globdata.colrprio:
            probspec = array([sum(globdata.maxmnumbpnts) / 2.])
            probsind = array([sum(globdata.maxmnumbpnts) / 2.])
        else:
            probspec = array([sum(globdata.maxmnumbpnts) / 2.])
            #probspec = array([sum(globdata.maxmnumbpnts) / 2.] * globdata.numbener)
            probsind = array([0.])
            
        globdata.probprop = concatenate((probfdfnnorm, probfdfnslop, probpsfipara, probnormback, probbrth, probdeth, probsplt, probmerg, problgal, probbgal, probspec, probsind))
        globdata.probprop /= sum(globdata.probprop)
       

def retr_pixlcnvt(globdata):
 
    path = os.environ["PCAT_DATA_PATH"] + '/pixlcnvt_%03d.p' % globdata.maxmgang
    if os.path.isfile(path):
        fobj = open(path, 'rb')
        globdata.pixlcnvt = cPickle.load(fobj)
        fobj.close()
    else:
        globdata.pixlcnvt = zeros(globdata.numbpixlheal, dtype=int)
        for k in range(globdata.indxpixlrofimarg.size):
            dist = retr_dist(globdata, globdata.lgalheal[globdata.indxpixlrofimarg[k]], globdata.bgalheal[globdata.indxpixlrofimarg[k]], globdata.lgalgrid, globdata.bgalgrid)
            globdata.pixlcnvt[globdata.indxpixlrofimarg[k]] = argmin(dist)

        fobj = open(path, 'wb')
        cPickle.dump(globdata.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()


def retr_strgangl(globdata):

    if globdata.trueinfo:
        if globdata.datatype == 'mock':
            globdata.truelabl = 'Mock'
        if globdata.datatype == 'inpt':
            if globdata.exprtype == 'ferm':
                globdata.truelabl = '3FGL'
            if globdata.exprtype == 'sdss':
                globdata.truelabl = 'Hubble'

    if globdata.exprtype == 'ferm':
        globdata.strganglunit = '[deg]'
    if globdata.exprtype == 'sdss':
        globdata.strganglunit = '[arcsec]'

    if globdata.regitype == 'igal':
        globdata.longlabl = '$l$'
        globdata.latilabl = '$b$'
    if globdata.regitype == 'ngal':
        globdata.longlabl = r'$\nu$'
        globdata.latilabl = r'$\mu$'
    if globdata.exprtype == 'ferm':
        globdata.longlabl += r' [$^\circ$]'
        globdata.latilabl += r' [$^\circ$]'
    if globdata.exprtype == 'sdss':
        globdata.longlabl += ' [arcsec]'
        globdata.latilabl += ' [arcsec]'

    
def retr_indxcube(globdata):

    globdata.indxcubefull = meshgrid(globdata.indxener, globdata.indxpixl, globdata.indxevtt, indexing='ij')
    globdata.indxcubefilt = meshgrid(globdata.indxenerincl, globdata.indxpixl, globdata.indxevttincl, indexing='ij')
    if globdata.pixltype == 'heal':
        globdata.indxcubeheal = meshgrid(globdata.indxenerinclfull, globdata.indxpixlrofi, globdata.indxevttinclfull, indexing='ij')


def retr_expo(globdata):
 
    if globdata.strgexpo == 'unit':
        globdata.expo = ones((globdata.numbener, globdata.numbpixl, globdata.numbevtt))
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/' + globdata.strgexpo
        globdata.expo = pf.getdata(path)

        if globdata.pixltype == 'heal':
            globdata.expo = globdata.expo[globdata.indxcubeheal]
        else:
            globdata.expo = globdata.expo.reshape((globdata.expo.shape[0], -1, globdata.expo.shape[-1]))
            
        globdata.expo = globdata.expo[globdata.indxcubefilt]
    

def init(globdata):
    
    # the number of processes (each with globdata.numbswep samples)
    if globdata.numbproc == None:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu':
            globdata.numbproc = 10
        else:
            globdata.numbproc = 1
        
    if globdata.exprtype == 'ferm':
        if globdata.strgback == None:
            globdata.strgback = ['fermisotflux.fits', 'fermfdfmflux_ngal.fits']
        if globdata.lablback == None:
            globdata.lablback = [r'$\mathcal{I}$', r'$\mathcal{D}$']
        if globdata.nameback == None:
            globdata.nameback = ['normisot', 'normfdfm']

    globdata.strgfluxunit = retr_strgfluxunit(globdata)
        
    # number of bins
    globdata.numbspec = 10
    globdata.numbbins = 10

    globdata.minmnumbpnts = 1

    globdata.numbback = len(globdata.strgback)
    globdata.indxback = arange(globdata.numbback)
    
    if globdata.numbback == 1:
        globdata.probprop[5] = 0.
        
    # PSF class axis
    globdata.numbevtt = globdata.indxevttincl.size
    globdata.indxevtt = arange(globdata.numbevtt)
    
    # energy axis
    globdata.numbener = globdata.indxenerincl.size
    globdata.indxenerinclbins = empty(globdata.numbener+1, dtype=int)
    globdata.indxenerinclbins[0:-1] = globdata.indxenerincl
    globdata.indxenerinclbins[-1] = globdata.indxenerincl[-1] + 1
    globdata.binsener = array([0.1, 0.3, 1., 3., 10., 100.])[globdata.indxenerinclbins]
    globdata.diffener = (roll(globdata.binsener, -1) - globdata.binsener)[0:-1]

    globdata.meanener = sqrt(roll(globdata.binsener, -1) * globdata.binsener)[0:-1]
    globdata.indxener = arange(globdata.numbener, dtype=int)
    
    if globdata.colrprio: 
        # temp
        globdata.indxenerfdfn = array([globdata.numbener / 2])
        globdata.minmspec = globdata.minmspec * (globdata.meanener[globdata.indxenerfdfn] / globdata.meanener)**2
        globdata.maxmspec = globdata.maxmspec * (globdata.meanener[globdata.indxenerfdfn] / globdata.meanener)**2
        if globdata.datatype == 'mock':
            globdata.mockfdfnslop = tile(globdata.mockfdfnslop, (1, globdata.numbener))
    else:
        globdata.indxenerfdfn = globdata.indxener

    globdata.enerfdfn = globdata.meanener[globdata.indxenerfdfn]
        
    if globdata.colrprio:
        globdata.indxenerprio = globdata.indxenerfdfn
    else:
        globdata.indxenerprio = globdata.indxener

    globdata.indxspecpivt = globdata.numbspec / 2
    
    # temp
    if globdata.exprtype == 'sdss':
        globdata.diffener = ones(globdata.numbener)
        
    # angular globdata.deviation
    globdata.numbangl = 100
    if globdata.exprtype == 'sdss':
        globdata.maxmangldisp = deg2rad(5. / 3600.) # [rad]
    if globdata.exprtype == 'ferm':
        globdata.maxmangldisp = deg2rad(5.) # [rad]
    globdata.angldisp = linspace(0., globdata.maxmangldisp, globdata.numbangl) # [rad]
    maxmangl = deg2rad(3.5 * globdata.maxmgang) # [rad]
    angl = linspace(0., maxmangl, globdata.numbangl) # [rad]

    if globdata.trueinfo:
        if globdata.datatype == 'mock':
            globdata.truelabl = 'Mock'
        if globdata.datatype == 'inpt':
            globdata.truelabl = '3FGL'
            
    if globdata.exprtype == 'ferm':
        globdata.strganglunit = '[deg]'
    if globdata.exprtype == 'sdss':
        globdata.strganglunit = '[arcsec]'
        
    if globdata.exprtype == 'ferm':
        retr_fermpsfn(globdata)

    # energy bin string
    globdata.enerstrg, globdata.binsenerstrg = retr_enerstrg(globdata)
    
    # PSF class string
    globdata.evttstrg = []
    for m in globdata.indxevtt:
        globdata.evttstrg.append('PSF%d' % globdata.indxevttincl[m])
        

    globdata.spltjcbnfact = log(2.**(2 - globdata.numbener))
    
    # number of components
    globdata.numbcomp = 2 + globdata.numbener
    if globdata.colrprio:
        globdata.numbcomp += 1
    globdata.numbcompcolr = 4
    globdata.jcbnsplt = 2.**(2 - globdata.numbener)
    
    # convenience factors for CDF and ICDF transforms
    globdata.factfdfnnorm = log(globdata.maxmfdfnnorm / globdata.minmfdfnnorm)
    globdata.factfdfnslop = arctan(globdata.maxmfdfnslop) - arctan(globdata.minmfdfnslop)
    if globdata.colrprio:
        globdata.factsind = arctan(globdata.maxmsind) - arctan(globdata.minmsind)
    

    if globdata.datatype == 'mock':

        # mock PSF parameters
        if globdata.mockpsfntype == 'singgaus':
            globdata.numbmockformpara = 1
        elif globdata.mockpsfntype == 'singking':
            globdata.numbmockformpara = 2 
        elif globdata.mockpsfntype == 'doubgaus':
            globdata.numbmockformpara = 3
        elif globdata.mockpsfntype == 'gausking':
            globdata.numbmockformpara = 4
        elif globdata.mockpsfntype == 'doubking':
            globdata.numbmockformpara = 5

        globdata.numbmockpsfiparaevtt = globdata.numbener * globdata.numbmockformpara
        globdata.numbmockpsfipara = globdata.numbmockpsfiparaevtt * globdata.numbevtt
        globdata.indxmockpsfipara = arange(globdata.numbmockpsfipara)   


    if globdata.regitype == 'igal':
        globdata.longlabl = '$l$'
        globdata.latilabl = '$b$'
    else:
        globdata.longlabl = r'$\nu$'
        globdata.latilabl = r'$\mu$'
        
    if globdata.exprtype == 'ferm':
        globdata.longlabl += r' [$^\circ$]'
        globdata.latilabl += r' [$^\circ$]'
    if globdata.exprtype == 'sdss':
        globdata.longlabl += ' [arcsec]'
        globdata.latilabl += ' [arcsec]'
    
    # construct the PSF model
    retr_psfimodl(globdata)

    # proposals
    retr_propmodl(globdata)
    
    # factors in the prior expression
    globdata.priofactlgalbgal = 2. * log(1. / 2. / globdata.maxmgang)
    globdata.priofactfdfnslop = globdata.numbener * log(1. / (arctan(globdata.maxmfdfnslop) - arctan(globdata.minmfdfnslop)))
    globdata.priofactfdfnnorm = log(1. / (log(globdata.maxmfdfnnorm) - log(globdata.minmfdfnnorm)))

    # sample vector indices  
    globdata.indxsampnumbpnts = arange(globdata.numbpopl)
    globdata.indxsampfdfnnorm = arange(globdata.numbpopl) + amax(globdata.indxsampnumbpnts) + 1
    globdata.indxsampfdfnslop = arange(globdata.numbpopl * globdata.numbener).reshape((globdata.numbpopl, globdata.numbener)) + amax(globdata.indxsampfdfnnorm) + 1
    globdata.indxsamppsfipara = arange(globdata.numbpsfipara) + amax(globdata.indxsampfdfnslop) + 1
    globdata.indxsampnormback = arange(globdata.numbback * globdata.numbener).reshape((globdata.numbback, globdata.numbener)) + amax(globdata.indxsamppsfipara) + 1

    globdata.fluxpivt = sqrt(globdata.minmspec * globdata.maxmspec)
    
    globdata.maxmnumbcomp = globdata.maxmnumbpnts * globdata.numbcomp
    globdata.indxsampcompinit = amax(globdata.indxsampnormback) + 1
    
    # maximum number of parameters
    globdata.maxmsampsize = int(globdata.indxsampcompinit + globdata.maxmnumbcomp * globdata.numbpopl)
    
    if globdata.numbburn == None:
        globdata.numbburn = globdata.numbswep / 5
    if globdata.factthin == None:
        globdata.factthin = min(globdata.maxmsampsize * 5, globdata.numbswep / 2)


    # run tag
    globdata.rtag = retr_rtag(globdata, None)
    
    # plots
    if globdata.makeplot:
        if os.uname()[1] == 'fink1.rc.fas.harvard.edu' and getpass.getuser() == 'tansu':
            plotfold = '/n/pan/www/tansu/png/pcat/'
        else:
            plotfold = os.environ["PCAT_DATA_PATH"] + '/png/'
        globdata.plotpath = plotfold + globdata.strgtime + '_' + globdata.rtag + '/'
        cmnd = 'mkdir -p ' + globdata.plotpath
        os.system(cmnd)

    # number of samples to be saved
    globdata.numbsamp = (globdata.numbswep - globdata.numbburn) / globdata.factthin

    # rescale the positional update scale
    globdata.stdvlbhl /= 2. * globdata.maxmgang

    # determine proposal probabilities
    globdata.probpropminm = copy(globdata.probprop)
    globdata.probpropmaxm = copy(globdata.probprop)
    
    ## do not allow death or merge when the number of PS is at its minimum or 
    ## births or splits when the number of PS is at its maximum
    globdata.probpropminm[[globdata.indxpropdeth, globdata.indxpropmerg]] = 0.
    globdata.probpropmaxm[[globdata.indxpropbrth, globdata.indxpropsplt]] = 0.
    globdata.probprop /= sum(globdata.probprop)
    globdata.probpropmaxm /= sum(globdata.probpropmaxm)
    globdata.probpropminm /= sum(globdata.probpropminm)
    
    # population indices
    globdata.indxpopl = arange(globdata.numbpopl)

    globdata.maxmgangmarg = globdata.maxmgang + globdata.margsize
    
    globdata.minmlgalmarg = -globdata.maxmgangmarg
    globdata.maxmlgalmarg = globdata.maxmgangmarg
    globdata.minmbgalmarg = -globdata.maxmgangmarg
    globdata.maxmbgalmarg = globdata.maxmgangmarg
    globdata.minmlgal = -globdata.maxmgang
    globdata.maxmlgal = globdata.maxmgang
    globdata.minmbgal = -globdata.maxmgang
    globdata.maxmbgal = globdata.maxmgang
    

    # input data
    if globdata.datatype == 'inpt':
        
        path = os.environ["PCAT_DATA_PATH"] + '/' + globdata.strgexpr
        exprflux = pf.getdata(path)

        globdata.indxenerinclfull = arange(exprflux.shape[0])
        globdata.indxevttinclfull = arange(exprflux.shape[2])

        if globdata.pixltype == 'heal':
            globdata.numbpixlheal = exprflux.shape[1]
            globdata.numbsideheal = int(sqrt(globdata.numbpixlheal / 12))
        else:
            globdata.numbsidecart = exprflux.shape[1]
            exprflux = exprflux.reshape((exprflux.shape[0], globdata.numbsidecart**2, exprflux.shape[3]))
            
    else:
        
        if globdata.exprtype == 'ferm':
            globdata.indxenerinclfull = arange(5)
            globdata.indxevttinclfull = arange(4)
         
    if globdata.pixltype == 'heal':
        globdata.numbpixlheal = globdata.numbsideheal**2 * 12
        globdata.apix = 4. * pi / globdata.numbpixlheal
    if globdata.pixltype == 'cart':
        globdata.binslgalcart = linspace(globdata.minmlgal, globdata.maxmlgal, globdata.numbsidecart + 1)
        globdata.binsbgalcart = linspace(globdata.minmbgal, globdata.maxmbgal, globdata.numbsidecart + 1)
        globdata.lgalcart = (globdata.binslgalcart[0:-1] + globdata.binslgalcart[1:]) / 2.
        globdata.bgalcart = (globdata.binsbgalcart[0:-1] + globdata.binsbgalcart[1:]) / 2.
        globdata.apix = deg2rad(2. * globdata.maxmgang / globdata.numbsidecart)**2
        
    # temp
    globdata.tracsamp = False
    
    # center of the ROI
    if globdata.regitype == 'igal':
        globdata.cntrlghp, globdata.cntrbghp = 0., 0.
    else:
        globdata.cntrlghp, globdata.cntrbghp = 0., 90.
    
    # plot settings
    ## marker opacity
    globdata.mrkralph = 0.8
    ## marker size
    globdata.minmmrkrsize = 50
    globdata.maxmmrkrsize = 500
    ## ROI
    globdata.exttrofi = [globdata.minmlgal, globdata.maxmlgal, globdata.minmbgal, globdata.maxmbgal]
    if globdata.exprtype == 'sdss':
        globdata.exttrofi *= 3600.
        globdata.frambndr = globdata.maxmgang * 3600.
        globdata.frambndrmarg = globdata.maxmgangmarg * 3600.
    else:
        globdata.frambndr = globdata.maxmgang
        globdata.frambndrmarg = globdata.maxmgangmarg
     

    # FDM normalization prior limits
    globdata.factnormback = log(globdata.maxmnormback / globdata.minmnormback)

    # sky coordinates
    globdata.binslbhl = linspace(-globdata.maxmgang, globdata.maxmgang, globdata.numbbins + 1)

    globdata.binsspec = zeros((globdata.numbener, globdata.numbspec + 1))
    for i in globdata.indxener:
        globdata.binsspec[i, :] = logspace(log10(globdata.minmspec[i]), log10(globdata.maxmspec[i]), globdata.numbspec + 1)
    globdata.meanspec = sqrt(globdata.binsspec[:, 1:] * globdata.binsspec[:, 0:-1])
    globdata.diffspec = globdata.binsspec[:, 1:] - globdata.binsspec[:, 0:-1]

    if globdata.colrprio:
        globdata.binssindunit = arange(globdata.numbbins + 1.) / globdata.numbbins
        globdata.meansindunit = (globdata.binssindunit[:-1] + globdata.binssindunit[1:]) / 2.
        globdata.binssind = icdf_atan(globdata.binssindunit, globdata.minmsind, globdata.factsind)
        globdata.meansind = icdf_atan(globdata.meansindunit, globdata.minmsind, globdata.factsind)
        globdata.diffsind = globdata.binssind[1:] - globdata.binssind[:-1]

    globdata.numbspecprox = 10
    globdata.indxspecprox = arange(globdata.numbspecprox)
    globdata.meanspecprox = globdata.meanspec # globdata.binsspec[globdata.indxenerfdfn, globdata.numbspec / 2]
        
    if globdata.exprtype == 'ferm':
        globdata.maxmangleval = empty(globdata.numbspecprox)
        for h in globdata.indxspecprox:
            frac = 1e-2 * amax(globdata.minmspec) / amax(globdata.meanspecprox[:, h])
            globdata.maxmangleval[h] = rad2deg(amax(retr_psfnwdth(globdata, globdata.fermpsfn, frac)))
    if globdata.exprtype == 'sdss':
        globdata.maxmangleval = 10. / 3600.

    # pizelization
    if globdata.pixltype == 'heal':
        
        lgalheal, bgalheal, globdata.numbsideheal, globdata.numbpixlheal, globdata.apix = tdpy.util.retr_healgrid(globdata.numbsideheal)

        globdata.indxpixlrofi = where((abs(lgalheal) < globdata.maxmgang) & (abs(bgalheal) < globdata.maxmgang))[0]
        globdata.indxpixlrofimarg = where((abs(lgalheal) < globdata.maxmgangmarg + 300. / globdata.numbsideheal) & \
                (abs(bgalheal) < globdata.maxmgangmarg + 300. / globdata.numbsideheal))[0]
        
        globdata.lgalgrid = lgalheal[globdata.indxpixlrofi]
        globdata.bgalgrid = bgalheal[globdata.indxpixlrofi]
        
        path = os.environ["PCAT_DATA_PATH"] + '/pixlcnvt_%03d.p' % globdata.maxmgang
        if os.path.isfile(path):
            fobj = open(path, 'rb')
            globdata.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            globdata.pixlcnvt = zeros(globdata.numbpixlheal, dtype=int)
            for k in range(globdata.indxpixlrofimarg.size):
                dist = retr_dist(globdata, lgalheal[globdata.indxpixlrofimarg[k]], bgalheal[globdata.indxpixlrofimarg[k]], globdata.lgalgrid, globdata.bgalgrid)
                globdata.pixlcnvt[globdata.indxpixlrofimarg[k]] = argmin(dist)

            fobj = open(path, 'wb')
            cPickle.dump(globdata.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()


            
    else:
        isidecart = arange(globdata.numbsidecart)
        temp = meshgrid(isidecart, isidecart, indexing='ij')
        globdata.bgalgrid = globdata.bgalcart[temp[1].flatten()]
        globdata.lgalgrid = globdata.lgalcart[temp[0].flatten()]
        
    globdata.numbpixl = globdata.lgalgrid.size
    globdata.indxpixl = arange(globdata.numbpixl)
    globdata.indxcubefull = meshgrid(globdata.indxener, globdata.indxpixl, globdata.indxevtt, indexing='ij')
    
    globdata.indxcubefilt = meshgrid(globdata.indxenerincl, globdata.indxpixl, globdata.indxevttincl, indexing='ij')
    
    globdata.numbpixlsave = min(1000, globdata.numbpixl)
    globdata.indxpixlsave = choice(arange(globdata.numbpixlsave), size=globdata.numbpixlsave)


    if globdata.pixltype == 'heal':
        globdata.indxcubeheal = meshgrid(globdata.indxenerinclfull, globdata.indxpixlrofi, globdata.indxevttinclfull, indexing='ij')
        

    if globdata.datatype == 'inpt' and globdata.pixltype == 'heal':
        exprflux = exprflux[globdata.indxcubeheal]


    if globdata.datatype == 'inpt':
        exprflux = exprflux[globdata.indxcubefilt]
 
    if globdata.datatype == 'mock' or globdata.exprtype == 'ferm':
        globdata.mockpsfipara = globdata.fermpsfipara

    # exposure
    if globdata.strgexpo == 'unit':
        globdata.expo = ones((globdata.numbener, globdata.numbpixl, globdata.numbevtt))
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/' + globdata.strgexpo
        globdata.expo = pf.getdata(path)

        if globdata.pixltype == 'heal':
            globdata.expo = globdata.expo[globdata.indxcubeheal]
        else:
            globdata.expo = globdata.expo.reshape((globdata.expo.shape[0], -1, globdata.expo.shape[-1]))
            
        globdata.expo = globdata.expo[globdata.indxcubefilt]
    

    # backgrounds
    globdata.backflux = []
    globdata.backfluxmean = []
    for c, backfluxstrg in enumerate(globdata.strgback):
        path = os.environ["PCAT_DATA_PATH"] + '/' + backfluxstrg
        backfluxtemp = pf.getdata(path)
        if globdata.pixltype == 'heal':
            backfluxtemp = backfluxtemp[globdata.indxcubeheal]
        else:
            backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        backfluxtemp = backfluxtemp[globdata.indxcubefilt]
        globdata.backflux.append(backfluxtemp)
        globdata.backfluxmean.append(mean(sum(backfluxtemp * globdata.expo, 2) / sum(globdata.expo, 2), 1))
                
    # count axis
    globdata.expotemp = mean(globdata.expo, 1)
    globdata.minmcnts = globdata.minmspec * amin(globdata.expotemp, 1) * globdata.diffener
    globdata.maxmcnts = globdata.maxmspec * amax(globdata.expotemp, 1) * globdata.diffener
    globdata.binscnts = zeros((globdata.numbener, globdata.numbspec + 1))
    for i in globdata.indxener:
        globdata.binscnts[i, :] = logspace(log10(globdata.minmcnts[i]), log10(globdata.maxmcnts[i]), globdata.numbspec + 1) # [1]
        
    # get 3FGL catalog
    if globdata.exprtype == 'ferm':
        retr_fgl3(globdata)

    # get count data
    ## input data
    if globdata.datatype == 'inpt':
        globdata.datacnts = exprflux * globdata.expo * globdata.apix * globdata.diffener[:, None, None] # [1]

    ## mock data
    if globdata.datatype == 'mock':

        if globdata.mocknumbpnts == None:
            globdata.mocknumbpnts = empty(globdata.numbpopl)
            for l in globdata.indxpopl:
                globdata.mocknumbpnts[l] = random_integers(globdata.minmnumbpnts, globdata.maxmnumbpnts[l])
        
        if globdata.mockfdfnslop == None:
            temp = rand(globdata.numbener * globdata.numbpopl)
            globdata.mockfdfnslop = icdf_atan(temp, globdata.minmfdfnslop, globdata.factfdfnslop)
            
        if globdata.mockpsfipara == None: 
            globdata.mockpsfntype = psfntpye
            numbmockpsfipara = globdata.numbpsfipara
            globdata.mockpsfipara = empty(numbmockpsfipara)
            for k in arange(numbmockpsfipara):
                globdata.mockpsfipara[k] = icdf_psfipara(globdata, rand(), k)
   
        if globdata.mocknormback == None:
            for c in globdata.indxback:
                globdata.mocknormback[c, :] = icdf_logt(rand(globdata.numbener), globdata.minmnormback[c], globdata.factnormback[c])

        mockcnts = [[] for l in globdata.indxpopl]
        mocklgal = [[] for l in globdata.indxpopl]
        mockbgal = [[] for l in globdata.indxpopl]
        mockspec = [[] for l in globdata.indxpopl]
        if globdata.colrprio:
            mocksind = [[] for l in globdata.indxpopl]
        for l in globdata.indxpopl:
            mocklgal[l] = icdf_self(rand(globdata.mocknumbpnts[l]), -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg)
            mockbgal[l] = icdf_self(rand(globdata.mocknumbpnts[l]), -globdata.maxmgangmarg, 2. * globdata.maxmgangmarg) 
            mockspec[l] = empty((globdata.numbener, globdata.mocknumbpnts[l]))
            for i in globdata.indxenerfdfn:
                mockspec[l][i, :] = icdf_spec(globdata, rand(globdata.mocknumbpnts[l]), globdata.mockfdfnslop[l, i], globdata.minmspec[i], globdata.maxmspec[i])
            if globdata.colrprio:
                mocksind[l] = icdf_atan(rand(globdata.mocknumbpnts[l]), globdata.minmsind, globdata.factsind)
                mockspec[l] = retr_spec(globdata, mockspec[l][globdata.indxenerfdfn, :].flatten(), mocksind[l])
            indxpixltemp = retr_indxpixl(globdata, mockbgal[l], mocklgal[l])
            mockcnts[l] = mockspec[l][:, :, None] * globdata.expo[:, indxpixltemp, :] * globdata.diffener[:, None, None]
        mockpntsflux = retr_pntsflux(globdata, concatenate(mocklgal), concatenate(mockbgal), concatenate(mockspec, axis=1), globdata.mockpsfipara, globdata.mockpsfntype)
        mocktotlflux = retr_rofi_flux(globdata, globdata.mocknormback, mockpntsflux, globdata.indxcubefull)
        mocktotlcnts = mocktotlflux * globdata.expo * globdata.apix * globdata.diffener[:, None, None] # [1]

        globdata.datacnts = zeros((globdata.numbener, globdata.numbpixl, globdata.numbevtt))
        for i in globdata.indxener:
            for k in range(globdata.numbpixl):
                for m in globdata.indxevtt:
                    globdata.datacnts[i, k, m] = poisson(mocktotlcnts[i, k, m])
              
        if globdata.trueinfo:
            
            globdata.truelgal = []
            globdata.truebgal = []
            globdata.truespec = []
            if globdata.colrprio:
                globdata.truesind = []
            for i in globdata.indxpopl:
                globdata.truelgal.append(mocklgal[l])
                globdata.truebgal.append(mockbgal[l])
                globdata.truespec.append(mockspec[l])
                if globdata.colrprio:
                    globdata.truesind.append(mocksind[l])
                    
            globdata.indxtruepntstimevari = [array([])] * globdata.numbpopl
                    
            globdata.truenumbpnts = globdata.mocknumbpnts
            globdata.truefdfnslop = globdata.mockfdfnslop
            globdata.truenormback = globdata.mocknormback
            globdata.truecnts = mockcnts
            
            globdata.truespec = []
            for l in globdata.indxpopl:
                globdata.truespectemp = empty((3, globdata.numbener, globdata.mocknumbpnts[l]))
                globdata.truespectemp[:] = mockspec[l][None, :, :]
                globdata.truespec.append(globdata.truespectemp)
                
            globdata.truepsfipara = globdata.mockpsfipara


    ## Real data
    # true data
    if globdata.trueinfo:
        if globdata.datatype == 'inpt':
            globdata.truenumbpnts = None
            globdata.truefdfnslop = None
            globdata.truenormback = None
    
            if globdata.exprtype == 'ferm':
                globdata.truenumbpnts = array([globdata.fgl3numbpntsrofi], dtype=int)
                globdata.truelgal = [globdata.fgl3lgal[globdata.indxfgl3rofi]]
                globdata.truebgal = [globdata.fgl3bgal[globdata.indxfgl3rofi]]
                globdata.truespec = [globdata.fgl3spec[:, :, globdata.indxfgl3rofi]]
                if globdata.colrprio:
                    globdata.truesind = [globdata.fgl3sind[globdata.indxfgl3rofi]]
                indxpixltemp = retr_indxpixl(globdata, globdata.truebgal[0], globdata.truelgal[0])
                spec = globdata.fgl3spec[0, :, globdata.indxfgl3rofi]
                # temp
                spec = spec.T
                globdata.truecnts = [spec[:, :, None] * globdata.expo[:, indxpixltemp, :] * globdata.diffener[:, None, None]]
                globdata.indxtruepntstimevari = [globdata.indxfgl3timevarirofi]
                globdata.truepsfipara = globdata.fermpsfipara
                
        if globdata.datatype == 'mock':
            globdata.truepsfn = retr_psfn(globdata, globdata.truepsfipara, globdata.indxener, globdata.angldisp, globdata.mockpsfntype)
        else:
            if globdata.exprtype == 'ferm':
                globdata.truepsfn = globdata.fermpsfn
            if globdata.exprtype == 'sdss':
                globdata.truepsfn = sdsspsfn
                
        truefwhm = 2. * retr_psfnwdth(globdata, globdata.truepsfn, 0.5)
        
        truebackcnts = []
        globdata.truesigm = []
        for l in globdata.indxpopl:
            indxpixltemp = retr_indxpixl(globdata, globdata.truebgal[l], globdata.truelgal[l])
            truebackcntstemp = zeros((globdata.numbener, globdata.truenumbpnts[l], globdata.numbevtt))
            for c in globdata.indxback:
                truebackcntstemp += globdata.backflux[c][:, indxpixltemp, :] * globdata.expo[:, indxpixltemp, :] * \
                    globdata.diffener[:, None, None] * pi * truefwhm[:, None, :]**2 / 4.
            truebackcnts.append(truebackcntstemp)
            globdata.truesigm.append(globdata.truecnts[l] / sqrt(truebackcntstemp))
        
    # sanity checks
    if amax(abs(globdata.datacnts - globdata.datacnts.astype(int)) / globdata.datacnts) > 1e-3:
        print 'Fractional counts!'
        
    if amin(globdata.datacnts) < 0.:
        print 'Negative counts!'

    globdata.datafluxmean = mean(sum(globdata.datacnts, 2) / sum(globdata.expo, 2), 1) / globdata.apix / globdata.diffener
    globdata.datacntsmean = mean(sum(globdata.datacnts, 2), 1)
    globdata.datacntssatu = ceil((amax(sum(globdata.datacnts, 2), 1) - globdata.datacntsmean) * 0.05 + globdata.datacntsmean)
    globdata.resicntssatu = ceil(globdata.datacntssatu * 0.5)
    
    # auxiliary variables for plots
    if globdata.pixltype == 'heal':
        for i in globdata.indxener:
            for m in globdata.indxevtt:
                globdata.datacntscarttemp = tdpy.util.retr_cart(globdata.datacnts[i, :, m], globdata.indxpixlrofi, numbsideinpt=globdata.numbsideheal, \
                    minmlgal=globdata.minmlgal, maxmlgal=globdata.maxmlgal, minmbgal=globdata.minmbgal, maxmbgal=globdata.maxmbgal)
                if i == 0 and m == 0:
                    globdata.datacntscart = zeros((globdata.datacntscarttemp.shape[0], globdata.datacntscarttemp.shape[1], globdata.numbener, globdata.numbevtt))
                globdata.datacntscart[:, :, i, m] = globdata.datacntscarttemp
        
    else:
        globdata.datacntscart = globdata.datacnts.reshape((globdata.numbener, globdata.numbsidecart, globdata.numbsidecart, globdata.numbevtt))
        globdata.datacntscart = swapaxes(swapaxes(globdata.datacntscart, 0, 2), 0, 1)

    for i in globdata.indxener:
        indxdatacntscartsatu = where(globdata.datacntscart[:, :, i, :] > globdata.datacntssatu[i])
        globdata.datacntscart[indxdatacntscartsatu[0],                               indxdatacntscartsatu[1], i, indxdatacntscartsatu[2]] = globdata.datacntssatu[i]

    # make a look-up table of nearby pixels for each pixel
    path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%03d.p' % globdata.maxmgang
    if os.path.isfile(path):
        print 'Retrieving previously computed pixel look-up table...'
        fobj = open(path, 'rb')
        globdata.indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        print 'Computing the look-up table...'
        globdata.indxpixlprox = [[] for h in range(globdata.numbspecprox)]
        for j in globdata.indxpixl:
            dist = retr_dist(globdata, globdata.lgalgrid[j], globdata.bgalgrid[j], globdata.lgalgrid, globdata.bgalgrid)
            for h in range(globdata.numbspecprox):
                globdata.indxpixlproxtemp = where(dist < deg2rad(globdata.maxmangleval[h]))[0]
                globdata.indxpixlprox[h].append(globdata.indxpixlproxtemp)
        fobj = open(path, 'wb')
        cPickle.dump(globdata.indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()

