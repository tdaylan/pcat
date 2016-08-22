# common imports
from __init__ import *

class gdatstrt(object):
    
    def __init__(self):
        pass


def retr_psfnwdth(gdat, psfn, frac):
    
    fwhm = zeros((gdat.numbener, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            indxanglgood = argsort(psfn[i, :, m])
            intpfwhm = max(frac * amax(psfn[i, :, m]), amin(psfn[i, :, m]))
            if intpfwhm > amin(psfn[i, indxanglgood, m]) and intpfwhm < amax(psfn[i, indxanglgood, m]):
                fwhm[i, m] = interp1d(psfn[i, indxanglgood, m], gdat.binsangl[indxanglgood])(intpfwhm)
    return fwhm


def retr_spec(gdat, flux, sind):
    
    if isscalar(flux):
        flux = array([flux])

    if isscalar(sind):
        sind = array([sind])

    spec = flux[None, :] * (gdat.meanener[:, None] / gdat.enerfluxdist)**(-sind[None, :])
    
    return spec


def retr_indx(gdat, indxpntsfull):    

    indxsamplgal = []
    indxsampbgal = []
    indxsampspec = []
    indxsampsind = []
    indxsampcomp = []
    for l in gdat.indxpopl:
        indxsamplgaltemp = gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:l]) + array(indxpntsfull[l], dtype=int) * gdat.numbcomp
        indxsamplgal.append(indxsamplgaltemp)
        indxsampbgal.append(indxsamplgaltemp + 1)
        indxsampspec.append(repeat((indxsamplgaltemp + 2)[None, :], gdat.numbener, 0) +repeat(arange(gdat.numbener), len(indxpntsfull[l])).reshape(gdat.numbener, -1))
        indxsampsind.append(indxsamplgaltemp + 2 + gdat.numbener)
        indxsampcomp.append(repeat(indxsamplgaltemp, gdat.numbcomp) + tile(arange(gdat.numbcomp, dtype=int), len(indxpntsfull[l])))

    return indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp


def retr_pntsflux(gdat, lgal, bgal, spec, psfipara, psfntype):
    
    numbpnts = lgal.size
    pntsfluxsing = empty((numbpnts, gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for k in range(numbpnts):
        
        # calculate the distance to all pixels from each point source
        dist = retr_angldistunit(gdat, lgal[k], bgal[k], gdat.indxpixl)

        indx = argsort(dist)
        dist = dist[indx]
        indxpixltemp = gdat.indxpixl[indx]
        
        # evaluate the PSF
        psfn = retr_psfn(gdat, psfipara, gdat.indxener, dist, psfntype)
        for i in gdat.indxener:
            for m in gdat.indxevtt:
                pntsfluxsing[k, i, indxpixltemp, m] = spec[i, k] * psfn[i, :, m]
    
    # sum contributions from all PS
    pntsflux = sum(pntsfluxsing, 0) 

    return pntsflux


def retr_rofi_flux(gdat, normback, pntsflux, tempindx):
    
    modlflux = pntsflux[tempindx]
    for c in gdat.indxback:
        modlflux += normback[c, :, None, None] * gdat.backflux[c][tempindx]        
    
    return modlflux


def cdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):

    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (gdat.maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    fluxunit = norm / (1. - fluxdistsloplowr) * fluxdistbrek**fluxdistsloplowr * (flux**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
    indxflux = where(flux >= fluxdistbrek)[0]
    
    if indxflux.size > 0:
        temp = norm * fluxdistbrek**fluxdistsloplowr / (1. - fluxdistsloplowr) * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
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


def icdf_flux_brok(gdat, fluxunit, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr):
   
    norm = 1. / (fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr)) / (1. - fluxdistsloplowr) + \
                 fluxdistbrek**fluxdistslopuppr * (gdat.maxmflux**(1. - fluxdistslopuppr) - fluxdistbrek**(1. - fluxdistslopuppr)) / (1. - fluxdistslopuppr))
    fluxunitbrek = norm / (1. - fluxdistsloplowr) * fluxdistbrek**fluxdistsloplowr * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
    flux = (fluxunit * (1. - fluxdistsloplowr) / norm / fluxdistbrek**fluxdistsloplowr + gdat.minmflux**(1. - fluxdistsloplowr))**(1. / (1. - fluxdistsloplowr))
    indxfluxunit = where(fluxunit >= fluxunitbrek)[0]
    
    if indxfluxunit.size > 0:
        temp = norm * fluxdistbrek**fluxdistsloplowr / (1. - fluxdistsloplowr) * (fluxdistbrek**(1. - fluxdistsloplowr) - gdat.minmflux**(1. - fluxdistsloplowr))
        flux[indxfluxunit] = ((fluxunit[indxfluxunit] - temp) * (1. - fluxdistslopuppr) / norm / fluxdistbrek**fluxdistslopuppr + \
                                                                                            fluxdistbrek**(1. - fluxdistslopuppr))**(1. / (1. - fluxdistslopuppr))

    return flux


def cdfn_flux_powr(gdat, flux, fluxdistslop):
        
    fluxunit = (flux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop)) / (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop))
        
    return fluxunit


def icdf_flux_powr(gdat, fluxunit, fluxdistslop):

    flux = (fluxunit * (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop)) + gdat.minmflux**(1. - fluxdistslop))**(1. / (1. - fluxdistslop))
    
    return flux


def pdfn_flux_powr(gdat, flux, fluxdistslop):
  
    norm = (1. - fluxdistslop) / (gdat.maxmflux**(1. - fluxdistslop) - gdat.minmflux**(1. - fluxdistslop))
    
    pdfn = norm * flux**(-fluxdistslop)
    
    return pdfn


def icdf_self(paraunit, minmpara, factpara):
    para = factpara * paraunit + minmpara
    return para


def cdfn_self(para, minmpara, factpara):
    paraunit = (para - minmpara) / factpara
    return paraunit


def icdf_eerr(paraunit, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    
    cdfnnormpara = paraunit * cdfnnormdiff + cdfnnormminm
    tranpara = sp.special.erfinv(2. * cdfnnormpara - 1.) * sqrt(2)
    para = tranpara * stdvpara + meanpara
   
    return para


def cdfn_eerr(para, meanpara, stdvpara, cdfnnormminm, cdfnnormdiff):
    tranpara = (para - meanpara) / stdvpara
    cdfnnormpara = 0.5 * (sp.special.erf(tranpara / sqrt(2.)) + 1.)
    paraunit = (cdfnnormpara - cdfnnormminm) / cdfnnormdiff
    return paraunit


def icdf_logt(paraunit, minmpara, factpara):
    para = exp(paraunit * factpara) * minmpara
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


def icdf_psfipara(gdat, psfiparaunit, thisindxpsfipara):

    minmpsfipara = gdat.minmpsfipara[thisindxpsfipara]
    factpsfipara = gdat.factpsfipara[thisindxpsfipara]
    scalpsfipara = gdat.scalpsfipara[thisindxpsfipara]
        
    if scalpsfipara == 'self':
        psfipara = icdf_self(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfipara = icdf_logt(psfiparaunit, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfipara = icdf_atan(psfiparaunit, minmpsfipara, factpsfipara)

    return psfipara


def cdfn_psfipara(gdat, psfipara, thisindxpsfipara):
    
    minmpsfipara = gdat.minmpsfipara[thisindxpsfipara]
    factpsfipara = gdat.factpsfipara[thisindxpsfipara]
    scalpsfipara = gdat.scalpsfipara[thisindxpsfipara]

    if scalpsfipara == 'self':
        psfiparaunit = cdfn_self(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'logt':
        psfiparaunit = cdfn_logt(psfipara, minmpsfipara, factpsfipara)
    if scalpsfipara == 'atan':
        psfiparaunit = cdfn_atan(psfipara, minmpsfipara, factpsfipara)
    
    return psfiparaunit


def retr_thisindxprop(gdat, gdatmodi):
    
    # choose the population to be modified
    gdatmodi.indxpoplmodi = choice(gdat.indxpopl)
    
    numbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]
    if numbpnts == gdat.maxmnumbpnts[gdatmodi.indxpoplmodi]:
        gdatmodi.thisindxprop = choice(gdat.indxprop, p=gdat.probpropmaxm)
    elif numbpnts == gdat.minmnumbpnts:
        gdatmodi.thisindxprop = choice(gdat.indxprop, p=gdat.probpropminm)
    else:
        gdatmodi.thisindxprop = choice(gdat.indxprop, p=gdat.probprop)

    if gdat.verbtype > 1:
        print inspect.stack()[0][3]
        print 'indxpoplmodi'
        print gdatmodi.indxpoplmodi
        print 'thisnumbpnts'
        print numbpnts
        print 'maxmnumbpnts'
        print gdat.maxmnumbpnts
        print


def retr_indxpixl(gdat, bgal, lgal):

    if gdat.pixltype == 'heal':
        
        if False:
            temp = where(gdat.pixlcnvt > -1)[0]
            print 'hey'
            print 'bgal'
            print bgal
            print 'lgal'
            print lgal
            print 'temp'
            print amin(temp)
            print amax(temp)
            print 'ang2pix(gdat.numbsideheal, pi / 2. - bgal, lgal)'
            print ang2pix(gdat.numbsideheal, pi / 2. - bgal, lgal)
            print 'gdat.pixlcnvt'
            print gdat.pixlcnvt.size
            print amin(gdat.pixlcnvt)
            print amax(gdat.pixlcnvt)
            print 
            print
            print
    
        indxpixl = gdat.pixlcnvt[ang2pix(gdat.numbsideheal, pi / 2. - bgal, lgal)]
        if (indxpixl == -1).any():  
            print 'pixlcnvt went negative!'
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


def show_memo(objt, name):

    if isinstance(objt, list):
        for k in len(objt):
            size = sys.getsizeof(objt[k]) / 2.**20
    else:
        listsize = []
        listattr = []
        for attr, valu in objt.__dict__.iteritems():
            listsize.append(sys.getsizeof(valu) / 2.**20)
            listattr.append(attr)
        size = array(listsize)
        attr = array(listattr)
        sizetotl = sum(size) 
        print 'Memory budget: %s' % name
        print 'Total size: %.4g MB' % sizetotl
        
        # sort the sizes to get the largest tail
        indxsizetemp = argsort(size)[::-1]
        size = size[indxsizetemp]
        attr = attr[indxsizetemp]
        print 'Largest 5:'
        for k in range(5):
            print '%s: %.4g MB' % (attr[k], size[k])
        print 


def retr_llik(gdat, gdatmodi, init=False):

    if init:

        gdatmodi.thisllik = gdat.datacnts * log(gdatmodi.thismodlcnts) - gdatmodi.thismodlcnts
        
    elif gdatmodi.thisindxprop >= gdat.indxproppsfipara:

        # load convenience variables
        timebegn = time.time()
        
        ## save time by not reloading all PS into modi structure for PSF parameter updates
        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            lgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsamplgal)]
            bgal = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampbgal)]
            spec = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdat.indxenermodi, :]]
        if gdatmodi.thisindxprop >= gdat.indxpropbrth:
            lgal = gdatmodi.modilgal[:gdat.numbmodipnts]
            bgal = gdatmodi.modibgal[:gdat.numbmodipnts]
            spec = gdatmodi.modispec[meshgrid(gdat.indxenermodi, arange(gdat.numbmodipnts), indexing='ij')]
        
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 0] = timefinl - timebegn
        
        # determine pixels over which to evaluate the log-likelihood
        timebegn = time.time()
        
        if gdatmodi.thisindxprop == gdat.indxpropnormback:
            gdat.indxpixlmodi = gdat.indxpixl
        if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfipara:
            thisindxpixlprox = []
            for k in range(gdat.numbmodipnts):
                
                # temp -- this may not work for extreme color PS!
                # take the flux at the pivot energy
                if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                    fluxtemp = gdatmodi.thissampvarb[concatenate(gdatmodi.thisindxsampspec, axis=1)[gdat.indxenerfluxdist, k]]
                else:
                    fluxtemp = gdatmodi.modispec[gdat.indxenerfluxdist, k]

                # find the flux index
                indxfluxproxtemp = amin(where(gdat.binsfluxprox - fabs(fluxtemp) > 0.)[0]) - 1
                indxpixltemp = retr_indxpixl(gdat, bgal[k], lgal[k])

                if False:
                    print 'hey'
                    print 'llik'
                    print 'indxpixltemp'
                    print indxpixltemp
                    print 'bgal[k]'
                    print bgal[k]
                    print 'lgal[k]'
                    print lgal[k]

                thisindxpixlprox.append(indxpixlprox[indxfluxproxtemp][indxpixltemp])
            gdat.indxpixlmodi = unique(concatenate(thisindxpixlprox))
        

        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 1] = timefinl - timebegn

        # construct the mesh grid for likelihood evaluation
        timebegn = time.time()
        
        if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
            gdat.indxcubemodi = meshgrid(gdat.indxenermodi, gdat.indxpixlmodi, gdat.indxevtt, indexing='ij')

        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 2] = timefinl - timebegn

        # update the model point source flux map, if needed
        timebegn = time.time()
        
        if gdatmodi.thisindxprop == gdat.indxproppsfipara or gdatmodi.thisindxprop >= gdat.indxpropbrth:
            
            # temp
            #if gdatmodi.thisindxprop == gdat.indxproppsfipara:
            #    gdatmodi.nextpntsflux[gdat.indxcubemodi] = 0.
            #else:
            gdatmodi.nextpntsflux[gdat.indxcubemodi] = copy(gdatmodi.thispntsflux[gdat.indxcubemodi])
                
            # evaluate the PSF for each PS over the set of data pixels to be updated
            # temp
            if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                numbrept = 2
            else:
                numbrept = 1
            for n in range(numbrept):

                # grab the PSF
                if n == 0:
                    psfnintp = gdatmodi.thispsfnintp
                else:
                    psfnintp = gdatmodi.nextpsfnintp

                for k in range(gdat.numbmodipnts):
                    
                    # calculate the distance to the pixels to be updated
                    dist = retr_angldistunit(gdat, lgal[k], bgal[k], thisindxpixlprox[k])

                    # interpolate the PSF
                    psfn = psfnintp(dist)
    
                    # add the contribution of the PS to the the proposed flux map
                    for i in range(gdat.indxenermodi.size):
                        if gdatmodi.thisindxprop == gdat.indxproppsfipara:
                            if n == 0:
                                spectemp = -spec[i, k]
                            else:
                                spectemp = spec[i, k]
                        else:
                            spectemp = spec[i, k]
                        gdatmodi.nextpntsflux[gdat.indxenermodi[i], thisindxpixlprox[k], :] += spectemp * psfn[gdat.indxenermodi[i], :, :]


        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 3] = timefinl - timebegn

        # update the total model flux map
        timebegn = time.time()
       
        indxtemp = meshgrid(gdat.indxback, gdat.indxenermodi, indexing='ij')

        if gdatmodi.thisindxprop == gdat.indxpropnormback:
            normback = gdatmodi.nextsampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdatmodi.thispntsflux
        if gdatmodi.thisindxprop == gdat.indxproppsfipara or gdatmodi.thisindxprop >= gdat.indxpropbrth:
            normback = gdatmodi.thissampvarb[gdat.indxsampnormback[indxtemp]]
            pntsflux = gdatmodi.nextpntsflux
        gdatmodi.nextmodlflux[gdat.indxcubemodi] = retr_rofi_flux(gdat, normback, pntsflux, gdat.indxcubemodi)

        # temp
        if False:
            for i in range(gdat.indxenermodi.size):
                test = zeros(gdat.numbpixl)
                if gdat.strgprop[gdatmodi.thisindxprop] == 'normback':
                    gdat.indxenermodi = array([gdat.indxenermodi])
                
                test = zeros(gdat.numbpixl)
                test[gdat.indxpixlmodi] = 1.
                tdpy.util.plot_maps('/Users/tansu/Desktop/test/indxpixlmodi_%d.pdf' % gdatmodi.cntrswep, test, indxpixlrofi=gdat.indxpixlrofi, \
                                                                                    numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                
                test = zeros(gdat.numbpixl)
                test[gdat.indxpixlmodi] = gdatmodi.thispntsflux[gdat.indxenermodi[i], gdat.indxpixlmodi, 0]
                tdpy.util.plot_maps('/Users/tansu/Desktop/test/thispnts_%d.pdf' % gdatmodi.cntrswep, test, indxpixlrofi=gdat.indxpixlrofi, \
                                                                                    numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                
                test = zeros(gdat.numbpixl)
                test[gdat.indxpixlmodi] = gdatmodi.nextpntsflux[gdat.indxenermodi[i], gdat.indxpixlmodi, 0]
                tdpy.util.plot_maps('/Users/tansu/Desktop/test/nextpnts_%d.pdf' % gdatmodi.cntrswep, test, indxpixlrofi=gdat.indxpixlrofi, \
                                                                                    numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                
                test = zeros(gdat.numbpixl)
                test[gdat.indxpixlmodi] = gdatmodi.thismodlflux[gdat.indxenermodi[i], gdat.indxpixlmodi, 0]
                tdpy.util.plot_maps('/Users/tansu/Desktop/test/thismodl_%d.pdf' % gdatmodi.cntrswep, test, indxpixlrofi=gdat.indxpixlrofi, \
                                                                                    numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                
                test = zeros(gdat.numbpixl)
                test[gdat.indxpixlmodi] = gdatmodi.nextmodlflux[gdat.indxenermodi[i], gdat.indxpixlmodi, 0]
                tdpy.util.plot_maps('/Users/tansu/Desktop/test/nextmodl_%d.pdf' % gdatmodi.cntrswep, test, indxpixlrofi=gdat.indxpixlrofi, \
                                                                                    numbpixl=gdat.numbpixlfull, pixltype=gdat.pixltype, \
                                                                                    minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, \
                                                                                    minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 4] = timefinl - timebegn

        # calculate the total model count map
        timebegn = time.time()
        
        gdatmodi.nextmodlcnts[gdat.indxcubemodi] = gdatmodi.nextmodlflux[gdat.indxcubemodi] * gdat.expo[gdat.indxcubemodi] * \
            gdat.apix * gdat.diffener[gdat.indxenermodi, None, None] # [1]
        
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 5] = timefinl - timebegn

        # calculate the likelihood
        timebegn = time.time()
        
        gdatmodi.nextllik[gdat.indxcubemodi] = gdat.datacnts[gdat.indxcubemodi] * log(gdatmodi.nextmodlcnts[gdat.indxcubemodi]) - gdatmodi.nextmodlcnts[gdat.indxcubemodi]
        
        timefinl = time.time()
        gdatmodi.listchrollik[gdatmodi.cntrswep, 6] = timefinl - timebegn

        if not isfinite(gdatmodi.nextllik[gdat.indxcubemodi]).any():
            warnings.warn('Log-likelihood went NAN!')
            
        gdatmodi.deltllik = sum(gdatmodi.nextllik[gdat.indxcubemodi] - gdatmodi.thisllik[gdat.indxcubemodi])
    
    else:
        gdatmodi.deltllik = 0.
        
    
def retr_backfwhmcnts(gdat, normback, fwhm):

    backfwhmcnts = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
    for c in gdat.indxback:
        backfwhmcnts += normback[c, :, None, None] * gdat.backflux[c] * gdat.expo * gdat.diffener[:, None, None] * pi * fwhm[:, None, :]**2 / 4.

    return backfwhmcnts


def retr_sigm(gdat, cnts, backfwhmcnts):
    
    sigm = cnts / sum(mean(backfwhmcnts, 1), 1)[:, None]

    return sigm


def retr_lpri(gdat, gdatmodi, init=False):
        
    if init:
            
        for l in gdat.indxpopl:
            thismeanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]]
            thisnumbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[l]]
            thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][gdat.indxenerfluxdist, :]]
            if gdat.bindprio:
        
                meanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[l]]
                if gdat.fluxdisttype[l] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[l]]
                    fluxhistmodl = meanpnts * pdfn_flux_powr(gdat, gdat.meanflux, fluxdistslop) * gdat.diffflux
                if gdat.fluxdisttype[l] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[l]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[l]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[l]]
                    fluxhistmodl = meanpnts * pdfn_flux_brok(gdat, gdat.meanflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr) * gdat.diffflux 
                
                fluxhist = histogram(thisflux, gdat.binsflux)[0]
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)

                gdatmodi.thislpri[l, :] = lprbpois
               
            else:
                        
                gdatmodi.thislpri[l, 0] = thisnumbpnts * log(thismeanpnts) - thismeanpnts - sp.special.gammaln(thisnumbpnts + 1)
               
                if gdat.fluxdisttype[l] == 'powr':
                    thisfluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[l]]
                    gdatmodi.thislpri[l, 1] = sum(log(pdfn_flux_powr(gdat, thisflux, thisfluxdistslop)))
                if gdat.fluxdisttype[l] == 'brok':
                    thisfluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[l]]
                    thisfluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[l]]
                    thisfluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[l]]
                    gdatmodi.thislpri[l, 1] = sum(log(pdfn_flux_brok(gdat, thisflux, thisfluxdistbrek, thisfluxdistsloplowr, thisfluxdistslopuppr)))
    
                if False:
                    print 'sort(thisflux)'
                    print sort(thisflux)
                    print 'log(pdfn_flux_powr(gdat, thisflux, thisfluxdistslop))'
                    print log(pdfn_flux_powr(gdat, sort(thisflux), thisfluxdistslop))
                    print 'gdatmodi.thislpri'
                    print gdatmodi.thislpri
                    print

            gdatmodi.nextlpri = copy(gdatmodi.thislpri)
                
    else:
        thismeanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
        thisnumbpnts = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]
        thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, :]]
        
        # determine if the log-prior needs to be recalculated
        boolupdtmeanpnts = gdatmodi.thisindxprop == gdat.indxpropmeanpnts
        boolupdtfluxdist = gdatmodi.thisindxprop == gdat.indxpropfluxdistslop or gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or \
                                             gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr
        boolupdtnumbpnts = gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg

        # compute the delta log-prior for the associated update
        if boolupdtmeanpnts or boolupdtfluxdist or boolupdtnumbpnts:

            # change number of PS
            if gdat.bindprio:
                
                meanpnts = gdatmodi.thissampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
                fluxhist = histogram(gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, :]], gdat.binsflux)[0] 
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
            
                
                # change mean number of PS
                if boolupdtmeanpnts:
                    meanpnts = gdatmodi.nextsampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
                
                # change flux distribution
                if boolupdtfluxdist:
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                        fluxdistslop = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                            fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                            fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        if  gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
                            fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                    if boolupdtnumbpnts:
                
                        if gdatmodi.thisindxprop == gdat.indxpropbrth:
                            fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfluxdist, 0], gdat.binsflux)[0]
                        elif gdatmodi.thisindxprop == gdat.indxpropdeth:
                            fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfluxdist, 0], gdat.binsflux)[0]
                        elif gdatmodi.thisindxprop == gdat.indxpropsplt:
                            fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfluxdist, 0], gdat.binsflux)[0]
                            fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfluxdist, 1:3], gdat.binsflux)[0]
                        elif gdatmodi.thisindxprop == gdat.indxpropmerg:
                            fluxhist -= histogram(-gdatmodi.modispec[gdat.indxenerfluxdist, 0:2], gdat.binsflux)[0]
                            fluxhist += histogram(gdatmodi.modispec[gdat.indxenerfluxdist, 2], gdat.binsflux)[0]
    
                # flux prior
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxhistmodl = meanpnts * pdfn_flux_powr(gdat, gdat.meanflux, fluxdistslop) * gdat.diffflux
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxhistmodl = meanpnts * pdfn_flux_brok(gdat, gdat.meanflux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr) * gdat.diffflux
            
                lprbpois = fluxhist * log(fluxhistmodl) - fluxhistmodl - sp.special.gammaln(fluxhist + 1)
                gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :] = lprbpois
            
                gdatmodi.deltlpri = sum(gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, :])

            else:
    
                if boolupdtnumbpnts or boolupdtmeanpnts:
                    if boolupdtnumbpnts:
                        nextnumbpnts = gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] = nextnumbpnts * log(thismeanpnts) - thismeanpnts - sp.special.gammaln(nextnumbpnts + 1)
                    else:
                        nextmeanpnts = gdatmodi.nextsampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] = thisnumbpnts * log(nextmeanpnts) - nextmeanpnts - sp.special.gammaln(thisnumbpnts + 1)
                    gdatmodi.deltlpri = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 0] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, 0]
                    
                if boolupdtfluxdist:
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                        nextfluxdistslop = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] = sum(log(pdfn_flux_powr(gdat, thisflux, nextfluxdistslop)))
                    if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                        nextfluxdistbrek = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                        nextfluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                        nextfluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                        gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] = sum(log(pdfn_flux_brok(gdat, thisflux, nextfluxdistbrek, nextfluxdistsloplowr, nextfluxdistslopuppr)))
                    gdatmodi.deltlpri = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, 1] - gdatmodi.thislpri[gdatmodi.indxpoplmodi, 1]
       
                    if False:
                        print 'sort(thisflux)'
                        print sort(thisflux)
                        print 'log(pdfn_flux_powr(gdat, thisflux, nextfluxdistslop))'
                        print log(pdfn_flux_powr(gdat, sort(thisflux), nextfluxdistslop))
                        print 'gdatmodi.thislpri'
                        print gdatmodi.thislpri
                        print

        else:
            gdatmodi.deltlpri = 0.
       
        
def retr_sampvarb(gdat, indxpntsfull, samp):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
    sampvarb = zeros_like(samp)
    sampvarb[gdat.indxsampnumbpnts] = samp[gdat.indxsampnumbpnts]
    for l in gdat.indxpopl:
        sampvarb[gdat.indxsampmeanpnts[l]] = icdf_logt(samp[gdat.indxsampmeanpnts[l]], gdat.minmmeanpnts[l], gdat.factmeanpnts[l])
        if gdat.fluxdisttype[l] == 'powr':
            sampvarb[gdat.indxsampfluxdistslop[l]] = icdf_atan(samp[gdat.indxsampfluxdistslop[l]], gdat.minmfluxdistslop[l], gdat.factfluxdistslop[l])
        if gdat.fluxdisttype[l] == 'brok':
            sampvarb[gdat.indxsampfluxdistbrek[l]] = icdf_logt(samp[gdat.indxsampfluxdistbrek[l]], gdat.minmfluxdistbrek[l], gdat.factfluxdistbrek[l])
            sampvarb[gdat.indxsampfluxdistsloplowr[l]] = icdf_atan(samp[gdat.indxsampfluxdistsloplowr[l]], gdat.minmfluxdistsloplowr[l], gdat.factfluxdistsloplowr[l])
            sampvarb[gdat.indxsampfluxdistslopuppr[l]] = icdf_atan(samp[gdat.indxsampfluxdistslopuppr[l]], gdat.minmfluxdistslopuppr[l], gdat.factfluxdistslopuppr[l])

    for k in gdat.indxpsfipara:
        sampvarb[gdat.indxsamppsfipara[k]] = icdf_psfipara(gdat, samp[gdat.indxsamppsfipara[k]], k)
    for c in gdat.indxback:
        sampvarb[gdat.indxsampnormback[c, :]] = icdf_logt(samp[gdat.indxsampnormback[c, :]], gdat.minmnormback[c], gdat.factnormback[c])
    
    for l in gdat.indxpopl:
        sampvarb[indxsamplgal[l]] = icdf_self(samp[indxsamplgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        sampvarb[indxsampbgal[l]] = icdf_self(samp[indxsampbgal[l]], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
        if gdat.fluxdisttype[l] == 'powr':
            sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_powr(gdat, samp[indxsampspec[l][gdat.indxenerfluxdist, :]], sampvarb[gdat.indxsampfluxdistslop[l]])
        if gdat.fluxdisttype[l] == 'brok':
            fluxunit = samp[indxsampspec[l][gdat.indxenerfluxdist[0], :]]
            fluxdistbrek = sampvarb[gdat.indxsampfluxdistbrek[l]]
            fluxdistsloplowr = sampvarb[gdat.indxsampfluxdistsloplowr[l]]
            fluxdistslopuppr = sampvarb[gdat.indxsampfluxdistslopuppr[l]]
            sampvarb[indxsampspec[l][gdat.indxenerfluxdist, :]] = icdf_flux_brok(gdat, fluxunit, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        sampvarb[indxsampsind[l]] = icdf_eerr(samp[indxsampsind[l]], gdat.sinddistmean[l], gdat.sinddiststdv[l], gdat.sindcdfnnormminm[l], gdat.sindcdfnnormdiff[l])
        sampvarb[indxsampspec[l]] = retr_spec(gdat, sampvarb[indxsampspec[l][gdat.indxenerfluxdist[0], :]], sampvarb[indxsampsind[l]])
    
    return sampvarb
    

def retr_maps(gdat, indxpntsfull, sampvarb):
    
    indxsamplgal, indxsampbgal, indxsampspec, indxsampsind, indxsampcomp = retr_indx(gdat, indxpntsfull)    
     
    listspectemp = []
    for l in gdat.indxpopl:
        listspectemp.append(sampvarb[indxsampspec[l]])

    pntsflux = retr_pntsflux(gdat, sampvarb[concatenate(indxsamplgal)], sampvarb[concatenate(indxsampbgal)], \
        concatenate(listspectemp, axis=1), sampvarb[gdat.indxsamppsfipara], gdat.psfntype)
    totlflux = retr_rofi_flux(gdat, sampvarb[gdat.indxsampnormback], pntsflux, gdat.indxcube)
    
    pntscnts = pntsflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    totlcnts = totlflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None]
    
    return pntsflux, pntscnts, totlflux, totlcnts


def retr_mrkrsize(gdat, flux):

    mrkrsize = (log(flux) - log(gdat.minmflux)) / (log(gdat.maxmflux) - log(gdat.minmflux)) * (gdat.maxmmrkrsize - gdat.minmmrkrsize) + gdat.minmmrkrsize
        
    return mrkrsize


def retr_fermpsfn(gdat):
   
    if False:
        reco = 8
    else:
        reco = 7

    if reco == 8:
        path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P7REP_SOURCE_V15_back.fits'
    irfn = pf.getdata(path, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = sqrt(minmener * maxmener)

    numbfermscalpara = 3
    numbfermformpara = 5
    
    fermscal = zeros((gdat.numbevtt, numbfermscalpara))
    fermform = zeros((gdat.numbener, gdat.numbevtt, numbfermformpara))
    
    parastrg = ['ntail', 'score', 'gcore', 'stail', 'gtail']
    for m in gdat.indxevtt:
        if reco == 8:
            irfn = pf.getdata(path, 1 + 3 * gdat.indxevttincl[m])
            fermscal[m, :] = pf.getdata(path, 2 + 3 * gdat.indxevttincl[m])['PSFSCALE']
        else:
            if m == 1:
                path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P7REP_SOURCE_V15_front.fits'
            elif m == 0:
                path = os.environ["PCAT_DATA_PATH"] + '/irfn/psf_P7REP_SOURCE_V15_back.fits'
            else:
                continue
            irfn = pf.getdata(path, 1)
            fermscal[m, :] = pf.getdata(path, 2)['PSFSCALE']
        for k in range(numbfermformpara):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)
        
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 0] = 1. / (1. + fermform[i, m, 0] * fermform[i, m, 3]**2 / fermform[i, m, 1]**2)

    # store the fermi PSF parameters
    gdat.fermpsfipara = zeros((gdat.numbener * numbfermformpara * gdat.numbevtt))
    for m in gdat.indxevtt:
        for k in range(numbfermformpara):
            indxfermpsfiparatemp = m * numbfermformpara * gdat.numbener + gdat.indxener * numbfermformpara + k
            gdat.fermpsfipara[indxfermpsfiparatemp] = fermform[:, m, k]

    # calculate the scale factor
    factener = (10. * gdat.meanener[:, None])**fermscal[None, :, 2]
    gdat.fermscalfact = sqrt((fermscal[None, :, 0] * factener)**2 + fermscal[None, :, 1]**2)
    
    # evaluate the PSF
    gdat.fermpsfn = retr_psfn(gdat, gdat.fermpsfipara, gdat.indxener, gdat.binsangl, 'doubking')


def retr_sdsspsfn(gdat):
   
    numbpsfiparaevtt = gdat.numbener * 3
    frac = array([[0.5, 0.5, 0.5]]).T
    sigc = array([[0.5, 1., 3.]]).T
    sigt = array([[0.5, 1., 3.]]).T
    gdat.sdsspsfipara = empty(numbpsfiparaevtt)
    gdat.sdsspsfn = retr_doubgaus(gdat.binsangl[None, :, None], frac[:, None, :], sigc[:, None, :], sigt[:, None, :])


def updt_samp(gdat, gdatmodi):
    
    if gdatmodi.thisindxprop == gdat.indxpropmeanpnts:
        gdatmodi.thissampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]]
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]
    
    # determine if a hyperparameter is to be updated
    thisbool = False
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistslop or gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or \
                                                gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
        thisbool = True
    
    # update the hyperparameters
    if thisbool:
        
        ## update the sample vector
        flux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]]
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
            fluxdistslop = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
            gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
            fluxunit = cdfn_flux_powr(gdat, flux, fluxdistslop)
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
            if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
                fluxdistbrek = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
            elif gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
                fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                fluxdistsloplowr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
            else:
                fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                fluxdistslopuppr = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
            fluxunit = cdfn_flux_brok(gdat, flux, fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)

        ## update the unit sample vector -- this is unique for hyperparameter updates
        gdatmodi.drmcsamp[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, :], -1] = fluxunit
        
        # update the prior register
        gdatmodi.thislpri[gdatmodi.indxpoplmodi, :] = gdatmodi.nextlpri[gdatmodi.indxpoplmodi, :]

    # proposals that change the likelihood
    if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
        gdatmodi.thisllik[gdat.indxcubemodi] = gdatmodi.nextllik[gdat.indxcubemodi]
        gdatmodi.thismodlcnts[gdat.indxcubemodi] = gdatmodi.nextmodlcnts[gdat.indxcubemodi]
        
    # PSF
    if gdatmodi.thisindxprop == gdat.indxproppsfipara:
        # temp
        gdatmodi.thispsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
     
    # background normalization
    if gdatmodi.thisindxprop == gdat.indxpropnormback:
        gdatmodi.thissampvarb[gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]] = \
            gdatmodi.nextsampvarb[gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]]
        
    # proposals that change the PS flux map
    if gdatmodi.thisindxprop >= gdat.indxpropbrth or gdatmodi.thisindxprop == gdat.indxproppsfipara:
        gdatmodi.thispntsflux[gdat.indxcubemodi] = copy(gdatmodi.nextpntsflux[gdat.indxcubemodi])

    # transdimensinal updates
    if gdatmodi.thisindxprop >= gdat.indxpropbrth and gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]
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
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompsind]] = gdatmodi.modisind[0]
        
    ## death
    if gdatmodi.thisindxprop == gdat.indxpropdeth:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdat.killindxpnts)
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdat.killindxpnts)

    ## split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:

        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].append(gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0])
        del gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0]
        
        ### update the components
        #### first component
        gdatmodi.thissampvarb[gdat.indxsampfrst] = gdatmodi.modilgal[1]
        gdatmodi.thissampvarb[gdat.indxsampfrst+1] = gdatmodi.modibgal[1]
        gdatmodi.thissampvarb[gdat.indxsampfrst+2:gdat.indxsampfrst+2+gdat.numbener] = gdatmodi.modispec[:, 1]
        #### second component
        gdatmodi.thissampvarb[gdat.indxsampseco] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdat.indxsampseco+1] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdat.indxsampseco+2:gdat.indxsampseco+2+gdat.numbener] = gdatmodi.modispec[:, 2]
        
    ## merge
    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        ### update the PS index lists
        gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi].remove(gdat.mergindxseco)
        gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi].append(gdat.mergindxseco)

        ### update the component
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcomplgal]] = gdatmodi.modilgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompbgal]] = gdatmodi.modibgal[2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompspec]] = gdatmodi.modispec[:, 2]
        gdatmodi.thissampvarb[gdatmodi.indxsampmodi[gdat.indxcompsind]] = gdatmodi.modisind[2]
        
    ## PS parameter proposals
    if gdatmodi.thisindxprop >= gdat.indxproplgal:  
        if gdatmodi.thisindxprop == gdat.indxproplgal:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modilgal[1]
        elif gdatmodi.thisindxprop == gdat.indxpropbgal:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modibgal[1]
        else:
            gdatmodi.thissampvarb[gdatmodi.indxsampmodispec] = gdatmodi.modispec[:, 1]
            if gdatmodi.thisindxprop == gdat.indxpropsind:
                gdatmodi.thissampvarb[gdatmodi.indxsampmodi] = gdatmodi.modisind[1]

    # update the unit sample vector
    if gdatmodi.indxsampmodi.size > 0:
        gdatmodi.drmcsamp[gdatmodi.indxsampmodi, 0] = gdatmodi.drmcsamp[gdatmodi.indxsampmodi, 1]


def retr_listpair(gdat, lgal, bgal):
    
    pairlist = []
    for k in range(lgal.size):
        indxpnts = k + 1 + where((lgal[k+1:] < lgal[k] + gdat.radispmrlbhl) & \
            (lgal[k+1:] > lgal[k] - gdat.radispmrlbhl) & (bgal[k+1:] < bgal[k] + gdat.radispmrlbhl) & (bgal[k+1:] > bgal[k] - gdat.radispmrlbhl))[0]
        for l in range(indxpnts.size):
            pairlist.append([k, indxpnts[l]])
            
    return pairlist


def retr_fermdata(gdat):
    
    path = os.environ["PCAT_DATA_PATH"] + '/gll_psc_v16.fit'

    fgl3 = pf.getdata(path)
    
    fgl3spectemp = stack((fgl3['Flux100_300'], fgl3['Flux300_1000'], fgl3['Flux1000_3000'], fgl3['Flux3000_10000'], fgl3['Flux10000_100000']))
    fgl3spectemp = fgl3spectemp[gdat.indxenerincl, :] / gdat.diffener[:, None]
    
    # sort the catalog in decreasing flux
    indxfgl3sort = argsort(fgl3spectemp[gdat.indxenerfluxdist[0], :])[::-1]

    fgl3spectemp = fgl3spectemp[:, indxfgl3sort]

    fgl3specstdvtemp = stack((fgl3['Unc_Flux100_300'], fgl3['Unc_Flux300_1000'], fgl3['Unc_Flux1000_3000'], fgl3['Unc_Flux3000_10000'], fgl3['Unc_Flux10000_100000']))
    fgl3specstdvtemp = fgl3specstdvtemp[gdat.indxenerincl, :, :] / gdat.diffener[:, None, None]
    fgl3specstdvtemp = fgl3specstdvtemp[:, indxfgl3sort, :]

    fgl3numbpntsfull = fgl3['glon'].size
    
    fgl3lgal = fgl3['glon'][indxfgl3sort]
    fgl3lgal = ((fgl3lgal - 180.) % 360.) - 180.

    fgl3bgal = fgl3['glat'][indxfgl3sort]
            
    fgl3axisstdv = (fgl3['Conf_68_SemiMinor'][indxfgl3sort] + fgl3['Conf_68_SemiMajor'][indxfgl3sort]) * 0.5
    fgl3anglstdv = deg2rad(fgl3['Conf_68_PosAng'][indxfgl3sort]) # [rad]
    fgl3lgalstdv = fgl3axisstdv * abs(cos(fgl3anglstdv))
    fgl3bgalstdv = fgl3axisstdv * abs(sin(fgl3anglstdv))

    fgl3sind = fgl3['Spectral_Index'][indxfgl3sort]
    
    fgl3strg = fgl3['Source_Name'][indxfgl3sort]
    fgl3strgclss = fgl3['CLASS1'][indxfgl3sort]
    fgl3strgassc = fgl3['ASSOC1'][indxfgl3sort]
    
    fgl3spectype = fgl3['SpectrumType'][indxfgl3sort]
    fgl3scur = fgl3['beta'][indxfgl3sort]
    fgl3scut = fgl3['Cutoff'][indxfgl3sort] * 1e-3
    
    fgl3timevari = fgl3['Variability_Index'][indxfgl3sort]
    
    fgl3spec = zeros((3, gdat.numbener, fgl3numbpntsfull))
    fgl3spec[0, :, :] = fgl3spectemp
    fgl3spec[1, :, :] = fgl3spectemp - fgl3specstdvtemp[:, :, 0]
    fgl3spec[2, :, :] = fgl3spectemp + fgl3specstdvtemp[:, :, 1]
    
    # adjust 3FGL positions according to the ROI center
    if gdat.regitype == 'ngal':
        rttr = hp.rotator.Rotator(rot=[0., 90., 0.], deg=True)
        fgl3bgal, fgl3lgal = rttr(pi / 2. - fgl3bgal, fgl3lgal)
        fgl3bgal = pi / 2. - fgl3bgal

    # select the 3FGL point sources in the ROI
    indxfgl3rofi = arange(fgl3lgal.size, dtype=int)
    for i in gdat.indxener:
        indxfgl3rofi = intersect1d(where((fgl3spec[0, i, :] > gdat.minmspec[i]) & (fgl3spec[0, i, :] < gdat.maxmspec[i]))[0], indxfgl3rofi)
    indxfgl3rofi = intersect1d(where((fabs(fgl3lgal) < gdat.maxmgangmarg) & (fabs(fgl3bgal) < gdat.maxmgangmarg))[0], indxfgl3rofi)

    # time variability
    indxfgl3vari = where(fgl3timevari[indxfgl3rofi] > 72.44)[0]
    
    # number of 3FGL PS in the ROI
    fgl3numbpnts = indxfgl3rofi.size

    # compute the 3FGL counts
    fgl3cnts = empty((gdat.numbener, fgl3numbpnts, gdat.numbevtt))
    for i in gdat.indxener:
        for m in gdat.indxevtt:
            # temp -- zero exposure pixels will give zero counts
            indxpixltemp = retr_indxpixl(gdat, fgl3bgal[indxfgl3rofi], fgl3lgal[indxfgl3rofi])
            fgl3cnts[i, :, m] = fgl3spec[0, i, indxfgl3rofi] * gdat.expo[i, indxpixltemp, m] * gdat.diffener[i]

    gdat.exprnumbpnts = fgl3numbpnts
    gdat.exprlgal = fgl3lgal[indxfgl3rofi]
    gdat.exprbgal = fgl3bgal[indxfgl3rofi]
    gdat.exprspec = fgl3spec[:, :, indxfgl3rofi]
    gdat.exprcnts = fgl3cnts
    gdat.exprsind = fgl3sind[indxfgl3rofi]
    gdat.exprstrg = fgl3strg[indxfgl3rofi]
    gdat.exprstrgclss = fgl3strgclss[indxfgl3rofi]
    gdat.exprstrgassc = fgl3strgassc[indxfgl3rofi]
    gdat.indxexprvari = indxfgl3vari
    
        

def retr_rtag(gdat, indxprocwork):
    
    if indxprocwork == None:
        rtag = 'AA_%d_%d_%d_%d_%s_%s_%s' % (gdat.numbproc, gdat.numbswep, gdat.numbburn, gdat.factthin, \
            gdat.datatype, gdat.regitype, gdat.psfntype)
    else:
        rtag = '%02d_%d_%d_%d_%d_%s_%s_%s' % (indxprocwork, gdat.numbproc, gdat.numbswep, gdat.numbburn, gdat.factthin, \
            gdat.datatype, gdat.regitype, gdat.psfntype)
        
    return rtag


def retr_gaus(gdat, gdatmodi, indxsamp, stdv):
    
    if gdat.fracrand > 0.:
        if rand() < gdat.fracrand:
            gdatmodi.drmcsamp[indxsamp, 1] = rand()
        else:
            gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, 0] + normal(scale=stdv)
    else:
        gdatmodi.drmcsamp[indxsamp, 1] = gdatmodi.drmcsamp[indxsamp, 0] + normal(scale=stdv)

        
def retr_angldistunit(gdat, lgal1, bgal1, indxpixltemp, retranglcosi=False):
    
    xaxi, yaxi, zaxi = retr_unit(lgal1, bgal1)
    anglcosi = gdat.xaxigrid[indxpixltemp] * xaxi + gdat.yaxigrid[indxpixltemp] * yaxi + gdat.zaxigrid[indxpixltemp] * zaxi
    if retranglcosi:
        return anglcosi
    else:
        angldist = arccos(anglcosi)
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


def retr_enerstrg(gdat):
    
    binsenerstrg = []
    if gdat.exprtype == 'ferm':
        enerstrg = []
        for i in gdat.indxener:
            binsenerstrg.append('%.3g GeV - %.3g GeV' % (gdat.binsener[i], gdat.binsener[i+1]))
            enerstrg.append('%.3g GeV' % gdat.meanener[i])
            
    if gdat.exprtype == 'sdss':
        gdat.binsenerstrg = ['i-band', 'r-band', 'g-band']
        enerstrg = ['i-band', 'r-band', 'g-band']
    
    if gdat.exprtype == 'chan':
        enerstrg = []
        for i in gdat.indxener:
            binsenerstrg.append('%.3g KeV - %.3g KeV' % (gdat.binsener[i] * 1e3, gdat.binsener[i+1] * 1e3))
            enerstrg.append('%.3g KeV' % (gdat.meanener[i] * 1e3))
    
    return enerstrg, binsenerstrg


def retr_prop(gdat, gdatmodi):
 
    gdatmodi.thisindxsamplgal, gdatmodi.thisindxsampbgal,  gdatmodi.thisindxsampspec, \
            gdatmodi.thisindxsampsind, gdatmodi.thisindxsampcomp = retr_indx(gdat, gdatmodi.thisindxpntsfull)
    
    if gdat.verbtype > 1:
        print 'retr_prop(): '

        print 'drmcsamp'
        print gdatmodi.drmcsamp
        
        print 'thissampvarb: '
        for k in range(gdatmodi.thissampvarb.size):
            if k == gdat.indxsampcompinit:
                print
            if k > gdat.indxsampcompinit and (k - gdat.indxsampcompinit) % gdat.numbcomp == 0:
                print
            print gdatmodi.thissampvarb[k]
        print
            
        print 'thisindxpntsfull: ', gdatmodi.thisindxpntsfull
        print 'thisindxpntsempt: ', gdatmodi.thisindxpntsempt  
        print 'thisindxsamplgal: ', gdatmodi.thisindxsamplgal
        print 'thisindxsampbgal: ', gdatmodi.thisindxsampbgal
        print 'thisindxsampspec: '
        print gdatmodi.thisindxsampspec
        print 'thisindxsampsind: ', gdatmodi.thisindxsampsind
        print 'thisindxsampcomp: ', gdatmodi.thisindxsampcomp
        print
        
    # hyperparameter changes
    # mean number of point sources
    if gdatmodi.thisindxprop == gdat.indxpropmeanpnts:
        gdatmodi.indxsampmodi = gdat.indxsampmeanpnts[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvmeanpnts)
        gdatmodi.nextsampvarb[gdat.indxsampmeanpnts] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
            gdat.minmmeanpnts[gdatmodi.indxpoplmodi], gdat.factmeanpnts[gdatmodi.indxpoplmodi])
        
    # flux distribution function shape
    ## flux distribution power law slope
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistslop:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfluxdistslop)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistslop[gdatmodi.indxpoplmodi], gdat.factfluxdistslop[gdatmodi.indxpoplmodi])
        gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampvarbmodi]), gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]))
        if gdat.verbtype > 1:
            print 'indxsampvarbmodi'
            print gdatmodi.indxsampvarbmodi
            print 'thissampvarb[gdat.indxsampfluxdistslop]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistslop]
            print 'nextsampvarb[gdat.indxsampfluxdistslop]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistslop]
    
    ## FDF break flux
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvflux)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistbrek[gdatmodi.indxpoplmodi], gdat.factfluxdistbrek[gdatmodi.indxpoplmodi])
        if gdat.verbtype > 1:
            print 'thissampvarb[gdat.indxsampfluxdistbrek]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek]
            print 'nextsampvarb[gdat.indxsampfluxdistbrek]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistbrek]
    
    ## FDF lower power law slope
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfluxdistslop)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistsloplowr[gdatmodi.indxpoplmodi], gdat.factfluxdistsloplowr[gdatmodi.indxpoplmodi])
        if gdat.verbtype > 1:
            print 'thissampvarb[gdat.indxsampfluxdistsloplowr]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr]
            print 'nextsampvarb[gdat.indxsampfluxdistsloplowr]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistsloplowr]
    
    ## FDF upper power law slope
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
        gdatmodi.indxsampvarbmodi = gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampvarbmodi, gdat.stdvfluxdistslop)
        gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]] = icdf_atan(gdatmodi.drmcsamp[gdatmodi.indxsampvarbmodi, -1], \
            gdat.minmfluxdistslopuppr[gdatmodi.indxpoplmodi], gdat.factfluxdistslopuppr[gdatmodi.indxpoplmodi])
        if gdat.verbtype > 1:
            print 'thissampvarb[gdat.indxsampfluxdistslopuppr]'
            print gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr]
            print 'nextsampvarb[gdat.indxsampfluxdistslopuppr]'
            print gdatmodi.nextsampvarb[gdat.indxsampfluxdistslopuppr]
   
    if gdatmodi.thisindxprop == gdat.indxpropfluxdistbrek or gdatmodi.thisindxprop == gdat.indxpropfluxdistsloplowr or gdatmodi.thisindxprop == gdat.indxpropfluxdistslopuppr:
        gdatmodi.indxsampmodi = concatenate((array([gdatmodi.indxsampvarbmodi]), gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], :]))

    # PSF parameter change 
    if gdatmodi.thisindxprop == gdat.indxproppsfipara:
        
        # index of the PSF parameter to change
        gdat.indxpsfiparamodi = choice(gdat.indxpsfipara)

        # the energy bin of the PS flux map to be modified
        gdat.indxenermodi = array([(gdat.indxpsfiparamodi % gdat.numbpsfiparaevtt) // gdat.numbformpara])
        
        # sample index to be modified
        gdatmodi.indxsampmodi = gdat.indxsamppsfipara[gdat.indxpsfiparamodi]
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvpsfipara)
        gdatmodi.nextsampvarb[gdat.indxsamppsfipara] = copy(gdatmodi.thissampvarb[gdat.indxsamppsfipara])
        gdatmodi.nextsampvarb[gdat.indxsamppsfipara[gdat.indxpsfiparamodi]] = \
            icdf_psfipara(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.indxpsfiparamodi)
            
        gdat.numbmodipnts = int(sum(gdatmodi.thissampvarb[gdat.indxsampnumbpnts]))
                    
        # construct the proposed PSF
        gdatmodi.nextpsfn = retr_psfn(gdat, gdatmodi.nextsampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.psfntype)
        
        temp = retr_psfn(gdat, gdatmodi.thissampvarb[gdat.indxsamppsfipara], gdat.indxener, gdat.binsangl, gdat.psfntype)
        gdatmodi.nextpsfnintp = interp1d(gdat.binsangl, gdatmodi.nextpsfn, axis=1)
        
        if gdat.verbtype > 1:
           
            print 'indxpsfiparamodi: ', gdat.indxpsfiparamodi
            print 'indxenermodi: ', gdat.indxenermodi
            print 'indxsampmodi: ', gdatmodi.indxsampmodi
            print 'thispsfipara'
            print gdatmodi.thissampvarb[gdat.indxsamppsfipara]
            print 'nextpsfipara'
            print gdatmodi.nextsampvarb[gdat.indxsamppsfipara]
            print 'thissampvarb[indxsampmodi]: ', gdatmodi.thissampvarb[gdatmodi.indxsampmodi]
            print 'nextpsfipara[indxsampmodi]: ', gdatmodi.nextsampvarb[gdatmodi.indxsampmodi]
            print

    # background changes
    if gdatmodi.thisindxprop == gdat.indxpropnormback:

        ## determine the sample index to be changed
        gdat.indxenermodi = choice(gdat.indxener)
        gdat.indxbackmodi = choice(gdat.indxback)
        gdatmodi.indxsampmodi = gdat.indxsampnormback[gdat.indxbackmodi, gdat.indxenermodi]
        
        ## propose
        retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvback)

        ## transform back from the unit space
        gdatmodi.nextsampvarb[gdat.indxsampnormback] = copy(gdatmodi.thissampvarb[gdat.indxsampnormback])
        gdatmodi.nextsampvarb[gdatmodi.indxsampmodi] = icdf_logt(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], \
                                                                                    gdat.minmnormback[gdat.indxbackmodi], gdat.factnormback[gdat.indxbackmodi])

        if gdat.verbtype > 1:
            print 'indxenermodi: ', gdat.indxenermodi
            print 'indxbackmodi: ', gdat.indxbackmodi
            print 'indxsampmodi: ', gdatmodi.indxsampmodi
            print 'thissampvarb[gdat.indxsampnormback]: ', gdatmodi.thissampvarb[gdat.indxsampnormback]
            print 'nextsampvarb[gdat.indxsampnormback]: ', gdatmodi.nextsampvarb[gdat.indxsampnormback]
            print

    # birth
    if gdatmodi.thisindxprop == gdat.indxpropbrth:

        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] + 1
    
        # initial sample index to add the new PS
        indxsampbrth = int(gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:gdatmodi.indxpoplmodi]) + gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0] * gdat.numbcomp)
        
        # sample auxiliary variables
        numbauxipara = gdat.numbcompcolr
        gdatmodi.auxipara = rand(numbauxipara)

        gdatmodi.drmcsamp[indxsampbrth:indxsampbrth+2, -1] = gdatmodi.auxipara[0:2]
        gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompflux, -1] = gdatmodi.auxipara[-2]
        gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompsind, -1] = gdatmodi.auxipara[-1]

        # sample indices to be modified
        gdatmodi.indxsampmodi = arange(indxsampbrth, indxsampbrth + gdat.numbcomp, dtype=int)

        # modification catalog
        gdat.numbmodipnts = 1
        gdatmodi.modilgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcomplgal, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        gdatmodi.modibgal[0] = icdf_self(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompbgal, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
        fluxunit = gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompflux, -1]
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
            fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
            modiflux = icdf_flux_powr(gdat, fluxunit, fluxdistslop)
        if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
            fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
            fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
            fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
            modiflux = icdf_flux_brok(gdat, array([fluxunit]), fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
        gdatmodi.modisind[0] = icdf_eerr(gdatmodi.drmcsamp[indxsampbrth+gdat.indxcompsind, -1], gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                    gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
        gdatmodi.modispec[:, 0] = retr_spec(gdat, modiflux, gdatmodi.modisind[0]).flatten()
    
        if gdat.verbtype > 1:
            print 'auxipara: ', gdatmodi.auxipara
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
    # kill
    if gdatmodi.thisindxprop == gdat.indxpropdeth:
        
        # change the number of PS
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # occupied PS index to be killed
        killindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be killed
        gdat.killindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][killindxindxpnts]
        
        # sample indices to be modified 
        gdatmodi.indxsampmodi = array([])
            
        # modification catalog
        gdat.numbmodipnts = 1
        gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][killindxindxpnts]]
        gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][killindxindxpnts]]
        gdatmodi.modispec[:, 0] = -gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, killindxindxpnts]]

        if gdat.verbtype > 1:
            print 'killindxpnts: ', gdat.killindxpnts
            print 'killindxindxpnts: ', killindxindxpnts
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print
            
  
    # split
    if gdatmodi.thisindxprop == gdat.indxpropsplt:
        
        gdat.numbmodipnts = 3
        
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] + 1
        
        # determine which point source to split
        thisindxindxpnts = arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int) 
        spltindxindxpnts = choice(thisindxindxpnts)
        spltindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][spltindxindxpnts]
    
        # update the sample vector
        gdat.indxsampfrst = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][spltindxindxpnts] * gdat.numbcomp
        indxfinlfrst = gdat.indxsampfrst + gdat.numbcomp
        gdat.indxsampseco = gdat.indxsampcompinit + gdat.maxmnumbcomp * gdatmodi.indxpoplmodi + gdatmodi.thisindxpntsempt[gdatmodi.indxpoplmodi][0] * gdat.numbcomp
        indxfinlseco = gdat.indxsampseco + gdat.numbcomp
        
        # determine the modified sample vector indices
        gdatmodi.indxsampmodi = concatenate((arange(gdat.indxsampfrst, indxfinlfrst, dtype=int), arange(gdat.indxsampseco, indxfinlseco, dtype=int)))
        
        thislgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][spltindxindxpnts]]
        thisbgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][spltindxindxpnts]]
        thisspec = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, spltindxindxpnts]]
        thisflux = thisspec[gdat.indxenerfluxdist[0]]
        thissind = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][spltindxindxpnts]]
        
        if gdat.verbtype > 1:
            print 'spltindxindxpnts: ', spltindxindxpnts
            print 'spltindxpnts: ', spltindxpnts
            print 'indxsampfrst: ', gdat.indxsampfrst
            print 'indxfinlfrst: ', indxfinlfrst
            print 'indxsampseco: ', gdat.indxsampseco
            print 'indxfinlseco: ', indxfinlseco
            print 'thislgal: ', gdat.anglfact * thislgal
            print 'thisbgal: ', gdat.anglfact * thisbgal
            print 'thisspec: ', thisspec
            print 'thisflux: ', thisflux
            print 'thissind: ', thissind
            
        # determine the new components
        # temp -- only valid for power-law energy spectrum
        gdatmodi.auxipara = empty(gdat.numbcompcolr)
        gdatmodi.auxipara[0] = rand()
        gdatmodi.auxipara[1] = rand() * gdat.radispmrlbhl
        gdatmodi.auxipara[2] = rand() * 2. * pi
        gdatmodi.auxipara[3] = icdf_eerr(rand(), gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                    gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
        
        if gdat.verbtype > 1:
            print 'auxipara[0]: ', gdatmodi.auxipara[0]
            print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
            print 'auxipara[2]: ', gdatmodi.auxipara[2]
            print 'auxipara[3]: ', gdatmodi.auxipara[3]
            print
            
        nextfluxfrst = gdatmodi.auxipara[0] * thisflux
        nextlgalfrst = thislgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1] * cos(gdatmodi.auxipara[2])
        nextbgalfrst = thislgal + (1. - gdatmodi.auxipara[0]) * gdatmodi.auxipara[1] * sin(gdatmodi.auxipara[2])
        nextsindfrst = thissind
        
        nextfluxseco = (1. - gdatmodi.auxipara[0]) * thisflux
        nextlgalseco = thislgal - gdatmodi.auxipara[0] * gdatmodi.auxipara[1] * cos(gdatmodi.auxipara[2])
        nextbgalseco = thislgal - gdatmodi.auxipara[0] * gdatmodi.auxipara[1] * sin(gdatmodi.auxipara[2])
        nextsindseco = gdatmodi.auxipara[3]
        
        if gdat.verbtype > 1:
            print 'nextlgalfrst: ', gdat.anglfact * nextlgalfrst
            print 'nextlgalseco: ', gdat.anglfact * nextlgalseco
            print 'nextbgalfrst: ', gdat.anglfact * nextbgalfrst
            print 'nextbgalseco: ', gdat.anglfact * nextbgalseco
            print 'nextfluxfrst: ', nextfluxfrst
            print 'nextfluxseco: ', nextfluxseco

        if fabs(nextlgalfrst) > gdat.maxmgangmarg or fabs(nextlgalseco) > gdat.maxmgangmarg or fabs(nextbgalfrst) > gdat.maxmgangmarg or fabs(nextbgalseco) > gdat.maxmgangmarg or \
                                    nextfluxfrst < gdat.minmflux or nextfluxseco < gdat.minmflux:
            gdatmodi.boolreje = True
                
        if not gdatmodi.boolreje:

            lgal = concatenate((array([nextlgalfrst, nextlgalseco]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], thislgal)))
            bgal = concatenate((array([nextbgalfrst, nextbgalseco]), setdiff1d(gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]], thisbgal)))
            listpair = retr_listpair(gdat, lgal, bgal)

            ## first new component
            gdatmodi.drmcsamp[gdat.indxsampfrst+gdat.indxcomplgal, -1] = cdfn_self(nextlgalfrst, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampfrst+gdat.indxcompbgal, -1] = cdfn_self(nextbgalfrst, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampfrst+gdat.indxcompflux, -1] = cdfn_flux_powr(gdat, nextfluxfrst, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdat.indxsampfrst+gdat.indxcompsind, -1] = cdfn_eerr(nextsindfrst, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            nextspecfrst = retr_spec(gdat, nextfluxfrst, nextsindfrst)

            ## second new component
            gdatmodi.drmcsamp[gdat.indxsampseco+gdat.indxcomplgal, -1] = cdfn_self(nextlgalseco, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampseco+gdat.indxcompbgal, -1] = cdfn_self(nextbgalseco, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampseco+gdat.indxcompflux, -1] = cdfn_flux_powr(gdat, nextfluxseco, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdat.indxsampseco+gdat.indxcompsind, -1] = cdfn_eerr(nextsindseco, gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                                gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            nextspecseco = retr_spec(gdat, nextfluxseco, nextsindseco)

            ## component to be removed
            gdatmodi.modilgal[0] = thislgal
            gdatmodi.modibgal[0] = thisbgal
            gdatmodi.modispec[:, 0] = -thisspec.flatten()

            ## first component to be added
            gdatmodi.modilgal[1] = nextlgalfrst
            gdatmodi.modibgal[1] = nextbgalfrst
            gdatmodi.modispec[:, 1] = nextspecfrst.flatten()

            # second component to be added
            gdatmodi.modilgal[2] = nextlgalseco
            gdatmodi.modibgal[2] = nextbgalseco
            gdatmodi.modispec[:, 2] = nextspecseco.flatten()

            tempindxindxpnts = setdiff1d(thisindxindxpnts, spltindxindxpnts)
            lgal = concatenate((gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][tempindxindxpnts]], array([nextlgalfrst]), array([nextlgalseco]))) 
            bgal = concatenate((gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][tempindxindxpnts]], array([nextbgalfrst]), array([nextbgalseco]))) 
            bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]]
            listpair = retr_listpair(gdat, lgal, bgal)
            numbpair = len(listpair)
        
    if gdatmodi.thisindxprop == gdat.indxpropmerg:
        
        gdat.numbmodipnts = 3
        
        gdatmodi.nextsampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] = gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]] - 1

        # determine the first PS to merge
        #dir2 = array([gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]], gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]]])
            
        lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi]]
        bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi]]
        listpair = retr_listpair(gdat, lgal, bgal)
        numbpair = len(listpair)
        
        if gdat.verbtype > 1:
            print 'lgal'
            print lgal
            print 'bgal'
            print bgal
            print 'pairlist'
            print listpair
           
        if numbpair == 0:
            gdatmodi.boolreje = True
        else:
            indxpairtemp = choice(arange(numbpair))
            mergindxindxpntsfrst = listpair[indxpairtemp][0]
            mergindxindxpntsseco = listpair[indxpairtemp][1]
  
        if not gdatmodi.boolreje:

            # fisrt PS index to be merged
            mergindxfrst = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]
            mergindxsampinit0 = gdat.indxsampcompinit + mergindxfrst * gdat.numbcomp

            # second PS index to be merged
            gdat.mergindxseco = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][mergindxindxpntsseco]
            mergindxsampinit1 = gdat.indxsampcompinit + gdat.mergindxseco * gdat.numbcomp

            # determine the modified sample vector indices
            gdat.indxsampfrst = gdat.indxsampcompinit + gdat.numbcomp * mergindxfrst
            indxfinlfrst = gdat.indxsampfrst + gdat.numbcomp
            gdat.indxsampseco = gdat.indxsampcompinit + gdat.numbcomp * gdat.mergindxseco
            indxfinlseco = gdat.indxsampseco + gdat.numbcomp

            gdatmodi.indxsampmodi = arange(gdat.indxsampfrst, indxfinlfrst)

            # indices of the PS to be merges
            mergindxpnts = sort(array([mergindxfrst, gdat.mergindxseco], dtype=int))

            thislgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            thisbgalfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            thisspecfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpntsfrst]]
            thissindfrst = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpntsfrst]]
            thisfluxfrst = thisspecfrst[gdat.indxenerfluxdist[0]]

            thislgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            thisbgalseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            thisspecseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][:, mergindxindxpntsseco]]
            thissindseco = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpntsseco]]
            thisfluxseco = thisspecseco[gdat.indxenerfluxdist[0]]

            # auxiliary parameters
            auxifrac = thisfluxfrst / (thisfluxfrst + thisfluxseco) 
            auxiradi = sqrt((thislgalseco - thislgalfrst)**2 + (thisbgalseco - thisbgalfrst)**2)
            auxiangl = arctan((thisbgalseco - thisbgalfrst) / (thislgalseco - thislgalfrst))
            auxisind = thissindseco

            # temp
            gdatmodi.auxipara = zeros(gdat.numbcompcolr)
            gdatmodi.auxipara[0] = auxifrac
            gdatmodi.auxipara[1] = auxiradi
            gdatmodi.auxipara[2] = auxiangl
            gdatmodi.auxipara[3] = thissindseco
            
            # merged PS
            nextflux = thisfluxfrst + thisfluxseco
            nextlgal = thislgalfrst + (1. - auxifrac) * (thislgalseco - thislgalfrst)
            nextbgal = thisbgalfrst + (1. - auxifrac) * (thisbgalseco - thisbgalfrst)
            nextsind = thissindfrst
            nextspec = retr_spec(gdat, nextflux, nextsind)

            # determine the unit variables for the merged PS
            gdatmodi.drmcsamp[gdat.indxsampfrst, -1] = cdfn_self(nextlgal, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampfrst+1, -1] = cdfn_self(nextbgal, -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.drmcsamp[gdat.indxsampfrst+2, -1] = cdfn_flux_powr(gdat, nextflux, gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]])
            gdatmodi.drmcsamp[gdat.indxsampfrst+3, -1] = gdatmodi.drmcsamp[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][mergindxindxpntsfrst], -2]

            gdat.numbmodipnts = 3
            ## first component to be merged
            gdatmodi.modilgal[0] = thislgalfrst
            gdatmodi.modibgal[0] = thisbgalfrst
            gdatmodi.modispec[:, 0] = -thisspecfrst.flatten()

            ## first component to be merged
            gdatmodi.modilgal[1] = thislgalseco
            gdatmodi.modibgal[1] = thisbgalseco
            gdatmodi.modispec[:, 1] = -thisspecseco.flatten()

            ## component to be added
            gdatmodi.modilgal[2] = nextlgal
            gdatmodi.modibgal[2] = nextbgal
            gdatmodi.modispec[:, 2] = nextspec.flatten()

            if gdat.verbtype > 1:
                print 'mergindxfrst: ', mergindxfrst
                print 'mergindxindxpntsfrst: ', mergindxindxpntsfrst
                print 'mergindxseco: ', gdat.mergindxseco
                print 'mergindxindxpntsseco: ', mergindxindxpntsseco
                print 'indxsampfrst: ', gdat.indxsampfrst
                print 'indxfinlfrst: ', indxfinlfrst
                print 'indxsampseco: ', gdat.indxsampseco
                print 'indxfinlseco: ', indxfinlseco
                print 'thislgalfrst: ', gdat.anglfact * thislgalfrst
                print 'thisbgalfrst: ', gdat.anglfact * thisbgalfrst
                print 'thislgalseco: ', gdat.anglfact * thislgalseco
                print 'thisbgalseco: ', gdat.anglfact * thisbgalseco
                print 'thisfluxfrst: ', thisfluxfrst
                print 'thisfluxseco: ', thisfluxseco
                print 'nextlgal: ', gdat.anglfact * nextlgal
                print 'nextbgal: ', gdat.anglfact * nextbgal
                print 'nextspec: ', nextspec
                print 'nextflux: ', nextflux
                print 'nextsind: ', nextsind
                print 'auxipara[0]: ', gdatmodi.auxipara[0]
                print 'auxipara[1]: ', gdat.anglfact * gdatmodi.auxipara[1]
                print 'auxipara[2]: ', gdatmodi.auxipara[2]
                print 'auxipara[3]: ', gdatmodi.auxipara[3]
                print

    # component change
    if gdatmodi.thisindxprop >= gdat.indxproplgal:     
        
        gdat.indxenermodi = gdat.indxener
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdatmodi.thisindxprop == gdat.indxproplgal:
                gdat.indxcompmodi = gdat.indxcomplgal
            else:
                gdat.indxcompmodi = gdat.indxcompbgal
        else:
            if gdatmodi.thisindxprop == gdat.indxpropflux:
                gdat.indxcompmodi = gdat.indxcompflux
            elif gdatmodi.thisindxprop == gdat.indxpropsind:
                gdat.indxcompmodi = gdat.indxcompsind
            
        # occupied PS index to be modified
        modiindxindxpnts = choice(arange(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]], dtype=int))
        
        # PS index to be modified
        modiindxpnts = gdatmodi.thisindxpntsfull[gdatmodi.indxpoplmodi][modiindxindxpnts]
        
        # initial sample index of the PS to be modified
        gdatmodi.indxsampmodiinit = gdat.indxsampcompinit + sum(gdat.maxmnumbcomp[:gdatmodi.indxpoplmodi]) + modiindxpnts * gdat.numbcomp
        
        # sample index to be modified
        gdatmodi.indxsampmodi = gdatmodi.indxsampmodiinit + gdat.indxcompmodi
        gdatmodi.indxsampmodispec = gdatmodi.indxsampmodiinit + 2 + gdat.indxener
        
        # propose
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvlbhl) 
        if gdatmodi.thisindxprop == gdat.indxpropflux:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvflux)
        if gdatmodi.thisindxprop == gdat.indxpropsind:
            retr_gaus(gdat, gdatmodi, gdatmodi.indxsampmodi, gdat.stdvsind)

        thisflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist[0], modiindxindxpnts]]
        thissind = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][modiindxindxpnts]]
        thisspec = retr_spec(gdat, thisflux, thissind)
            
        gdat.numbmodipnts = 2
        gdatmodi.modispec[:, 0] = -thisspec.flatten()
        if gdatmodi.thisindxprop == gdat.indxproplgal or gdatmodi.thisindxprop == gdat.indxpropbgal:
            if gdat.indxcompmodi == 0:
                gdatmodi.modilgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
                gdatmodi.modilgal[1] = icdf_self(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
                gdatmodi.modibgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            else:
                gdatmodi.modilgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
                gdatmodi.modibgal[0] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
                gdatmodi.modibgal[1] = icdf_self(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
            gdatmodi.modispec[:, 1] = thisspec.flatten()
        else:
            gdatmodi.modilgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            gdatmodi.modibgal[:2] = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            if gdatmodi.thisindxprop == gdat.indxpropflux:
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'powr':
                    fluxdistslop = gdatmodi.thissampvarb[gdat.indxsampfluxdistslop[gdatmodi.indxpoplmodi]]
                    modiflux = icdf_flux_powr(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], fluxdistslop)
                if gdat.fluxdisttype[gdatmodi.indxpoplmodi] == 'brok':
                    fluxdistbrek = gdatmodi.thissampvarb[gdat.indxsampfluxdistbrek[gdatmodi.indxpoplmodi]]
                    fluxdistsloplowr = gdatmodi.thissampvarb[gdat.indxsampfluxdistsloplowr[gdatmodi.indxpoplmodi]]
                    fluxdistslopuppr = gdatmodi.thissampvarb[gdat.indxsampfluxdistslopuppr[gdatmodi.indxpoplmodi]]
                    modiflux = icdf_flux_brok(gdat, gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], fluxdistbrek, fluxdistsloplowr, fluxdistslopuppr)
                gdatmodi.modisind[1] = gdatmodi.thissampvarb[gdatmodi.thisindxsampsind[gdatmodi.indxpoplmodi][modiindxindxpnts]]
            else:
                modiflux = gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[gdatmodi.indxpoplmodi][gdat.indxenerfluxdist, modiindxindxpnts]]
                gdatmodi.modisind[1] = icdf_eerr(gdatmodi.drmcsamp[gdatmodi.indxsampmodi, -1], gdat.sinddistmean[gdatmodi.indxpoplmodi], \
                        gdat.sinddiststdv[gdatmodi.indxpoplmodi], gdat.sindcdfnnormminm[gdatmodi.indxpoplmodi], gdat.sindcdfnnormdiff[gdatmodi.indxpoplmodi])
            gdatmodi.modispec[:, 1] = retr_spec(gdat, modiflux, gdatmodi.modisind[1]).flatten()

        if gdat.verbtype > 1:
            print 'modilgal: ', gdatmodi.modilgal
            print 'modibgal: ', gdatmodi.modibgal
            print 'modispec: '
            print gdatmodi.modispec
            print 'indxcompmodi: ', gdat.indxcompmodi
            print 'modiindxindxpnts: ', modiindxindxpnts
            print 'modiindxpnts: ', modiindxpnts

    if gdat.verbtype > 1:
        print 'indxpoplmodi'
        print gdatmodi.indxpoplmodi
        if not (gdatmodi.boolreje or gdatmodi.thisindxprop == gdat.indxpropdeth):
            print 'indxsampmodi'
            print gdatmodi.indxsampmodi
            print 'drmcsamp[gdatmodi.indxsampmodi, :]'
            print gdatmodi.drmcsamp[gdatmodi.indxsampmodi, :]
        
    # energy bin in which to evaluate the log-likelihood
    if gdat.indxpropbrth <= gdatmodi.thisindxprop <= gdat.indxpropmerg:
        gdat.indxenermodi = gdat.indxener

    if gdat.verbtype > 1:
        if gdatmodi.thisindxprop >= gdat.indxproppsfipara:
            print 'indxenermodi: ', gdat.indxenermodi

    # auxiliary variable density fraction and jacobian
    if (gdatmodi.thisindxprop == gdat.indxpropsplt or gdatmodi.thisindxprop == gdat.indxpropmerg) and not gdatmodi.boolreje:

        combfact = log(gdatmodi.thissampvarb[gdat.indxsampnumbpnts[gdatmodi.indxpoplmodi]]**2 / numbpair)
        if gdatmodi.thisindxprop == gdat.indxpropsplt:
            thisjcbnfact = gdat.spltjcbnfact
            thiscombfact = combfact 
        else:
            thisjcbnfact = -gdat.spltjcbnfact
            thiscombfact = -combfact 

        gdatmodi.laccfrac = thisjcbnfact + thiscombfact
        gdatmodi.listnumbpair[gdatmodi.cntrswep] = numbpair
        gdatmodi.listjcbnfact[gdatmodi.cntrswep] = thisjcbnfact
        gdatmodi.listcombfact[gdatmodi.cntrswep] = thiscombfact
        gdatmodi.listauxipara[gdatmodi.cntrswep, :] = gdatmodi.auxipara
        gdatmodi.listlaccfrac[gdatmodi.cntrswep] = gdatmodi.laccfrac

        if gdat.verbtype > 1:
            print 'thisjcbnfact'
            print thisjcbnfact
            print 'thiscombfact'
            print thiscombfact
            print 'laccfrac'
            print gdatmodi.laccfrac
            print 'listpair'
            print listpair
            print

    else:
        gdatmodi.laccfrac = 0.  
   
    # temp
    # define the index to be changed in the sample vector if a hyperparameter is being updated
    if gdatmodi.thisindxprop != gdat.indxpropfluxdistslop and gdatmodi.thisindxprop != gdat.indxpropfluxdistbrek and \
                                            gdatmodi.thisindxprop != gdat.indxpropfluxdistsloplowr and gdatmodi.thisindxprop != gdat.indxpropfluxdistslopuppr:
        gdatmodi.indxsampvarbmodi = gdatmodi.indxsampmodi

        
def retr_psfn(gdat, psfipara, indxenertemp, thisangl, psfntype):

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
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        scalangl = thisangl[None, :, None]

    indxpsfiparatemp = numbformpara * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    
    if psfntype == 'singgaus':
        sigc = psfipara[indxpsfiparatemp]
        sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
    elif psfntype == 'singking':
        sigc = psfipara[indxpsfiparatemp]
        gamc = psfipara[indxpsfiparatemp+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        psfn = retr_singking(scalangl, sigc, gamc)
        
    elif psfntype == 'doubgaus':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        psfn = retr_doubgaus(scalangl, frac, sigc, sigt)

    elif psfntype == 'gausking':
        frac = psfipara[indxpsfiparatemp]
        sigc = psfipara[indxpsfiparatemp+1]
        sigt = psfipara[indxpsfiparatemp+2]
        gamt = psfipara[indxpsfiparatemp+3]
        frac = frac[:, None, :]
        sigc = sigc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        psfn = retr_gausking(scalangl, frac, sigc, sigt, gamt)
        
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
        psfn = retr_doubking(scalangl, frac, sigc, gamc, sigt, gamt)
    
    # normalize the PSF
    psfn /= 2. * pi * trapz(psfn * sin(thisangl[None, :, None]), thisangl, axis=1)[:, None, :]
    
    # temp
    if True and (gdat.strgcnfg == 'pcat_ferm_expr_ngal' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp1' or \
                                                            gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp2' or gdat.strgcnfg == 'pcat_ferm_expr_ngal_cmp3'):
        print 'CORRECTING THE PSF.'
        tempcorr = array([1., 0.8, 0.8])
        psfn *= tempcorr[:, None, None]

    return psfn


def retr_psfimodl(gdat):
    
    # PSF parameters
    if gdat.psfntype == 'singgaus':
        gdat.numbformpara = 1
    elif gdat.psfntype == 'singking':
        gdat.numbformpara = 2 
    elif gdat.psfntype == 'doubgaus':
        gdat.numbformpara = 3
    elif gdat.psfntype == 'gausking':
        gdat.numbformpara = 4
    elif gdat.psfntype == 'doubking':
        gdat.numbformpara = 5
       
    gdat.indxformpara = arange(gdat.numbformpara) 
    gdat.numbpsfiparaevtt = gdat.numbener * gdat.numbformpara
    gdat.numbpsfipara = gdat.numbpsfiparaevtt * gdat.numbevtt
    gdat.indxpsfipara = arange(gdat.numbpsfipara)

    minmformpara = zeros(gdat.numbformpara)
    maxmformpara = zeros(gdat.numbformpara)
    factformpara = zeros(gdat.numbformpara)
    scalformpara = zeros(gdat.numbformpara, dtype=object)
    scalformpara = zeros(gdat.numbformpara, dtype=object)
    if gdat.exprtype == 'ferm':
        minmanglpsfn = 0.1 # deg2rad(0.0001)
        maxmanglpsfn = 5. # deg2rad(5.)
        #minmanglpsfn = 0.01
        #maxmanglpsfn = 3.
        # temp
        minmgamm = 2.
        maxmgamm = 20.
    if gdat.exprtype == 'chan' or gdat.exprtype == 'sdss':
        minmanglpsfn = 0.5 / gdat.anglfact
        maxmanglpsfn = 2. / gdat.anglfact 
    minmpsfnfrac = 0.
    maxmpsfnfrac = 1.
    
    if gdat.psfntype == 'singgaus':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        scalformpara[0] = 'logt'
        strgformpara = [r'$\sigma']
    elif gdat.psfntype == 'singking':
        minmformpara[0] = minmanglpsfn
        maxmformpara[0] = maxmanglpsfn
        minmformpara[1] = minmgamm
        maxmformpara[1] = maxmgamm
        scalformpara[0] = 'logt'
        scalformpara[1] = 'atan'
        strgformpara = [r'$\sigma', r'$\gamma']
    elif gdat.psfntype == 'doubgaus':
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
    elif gdat.psfntype == 'gausking':
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
    elif gdat.psfntype == 'doubking':
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

    for k in gdat.indxformpara:
        if scalformpara[k] == 'self':
            factformpara[k] = maxmformpara[k] - minmformpara[k]
        if scalformpara[k] == 'logt':
            factformpara[k] = log(maxmformpara[k] / minmformpara[k])
        if scalformpara[k] == 'atan':
            factformpara[k] = arctan(maxmformpara[k]) - arctan(minmformpara[k])
            
    gdat.minmpsfipara = tile(tile(minmformpara, gdat.numbener), gdat.numbevtt)
    gdat.maxmpsfipara = tile(tile(maxmformpara, gdat.numbener), gdat.numbevtt)
    gdat.scalpsfipara = tile(tile(scalformpara, gdat.numbener), gdat.numbevtt)
    gdat.factpsfipara = tile(tile(factformpara, gdat.numbener), gdat.numbevtt)
   
    # PSF parameter strings
    gdat.strgpsfipara = [strgformpara[k] + '^{%d%d}$' % (gdat.indxenerincl[i], gdat.indxevttincl[m]) \
        for m in gdat.indxevtt for i in gdat.indxener for k in gdat.indxformpara]
    gdat.indxpsfiparainit = (gdat.indxevtt[:, None] * gdat.numbpsfiparaevtt + gdat.indxener[None, :] * gdat.numbformpara).flatten()
    for k in arange(gdat.indxpsfiparainit.size):
        if gdat.psfntype == 'singgaus' or gdat.psfntype == 'singking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]] += ' ' + gdat.strganglunit
        elif gdat.psfntype == 'doubgaus' or gdat.psfntype == 'gausking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+1] += ' ' + gdat.strganglunit
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+2] += ' ' + gdat.strganglunit
        elif gdat.psfntype == 'doubking':
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+1] += ' ' + gdat.strganglunit
            gdat.strgpsfipara[gdat.indxpsfiparainit[k]+3] += ' ' + gdat.strganglunit
    gdat.indxpsfipara = arange(gdat.numbpsfipara)


def retr_unit(lgal, bgal):

    xaxi = cos(bgal) * cos(lgal)
    yaxi = -cos(bgal) * sin(lgal)
    zaxi = sin(bgal)

    return xaxi, yaxi, zaxi


def retr_propmodl(gdat):

    gdat.strgprop = []
    cntr = tdpy.util.cntr()
    
    ## mean number of PS
    gdat.indxpropmeanpnts = cntr.incr()
    gdat.strgprop.append('meanpnts')

    ## flux distribution
    gdat.indxpropfluxdistslop = cntr.incr()
    gdat.strgprop.append('fluxdistslop')
    gdat.indxpropfluxdistbrek = cntr.incr()
    gdat.strgprop.append('fluxdistbrek')
    gdat.indxpropfluxdistsloplowr = cntr.incr()
    gdat.strgprop.append('fluxdistsloplowr')
    gdat.indxpropfluxdistslopuppr = cntr.incr()
    gdat.strgprop.append('fluxdistslopuppr')

    # PSF parameters
    gdat.indxproppsfipara = cntr.incr()
    gdat.strgprop.append('psfipara')
    
    # background normalization
    gdat.indxpropnormback = cntr.incr()
    gdat.strgprop.append('normback')
    
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
    
    # sind
    gdat.strgprop.append('sind')
    gdat.indxpropsind = cntr.incr()

    gdat.numbprop = len(gdat.strgprop)
    gdat.indxprop = arange(gdat.numbprop)

    if gdat.probprop == None:
            
        if gdat.boolpropfluxdist:
            probmeanpnts = array([1.])
            probfluxdistslop = array([1.])
            probfluxdistbrek = array([1.])
            probfluxdistsloplowr = array([1.])
            probfluxdistslopuppr = array([1.])
        else:
            probmeanpnts = array([0.])
            probfluxdistslop = array([0.])
            probfluxdistbrek = array([0.])
            probfluxdistsloplowr = array([0.])
            probfluxdistslopuppr = array([0.])

        if gdat.boolproppsfn:
            probpsfipara = array([1.]) * gdat.numbpsfipara
        else:
            probpsfipara = array([0.])
        probnormback = array([1.])
        
        probbrth = array([0.2 * sum(gdat.maxmnumbpnts) / 2.])
        probdeth = array([0.2 * sum(gdat.maxmnumbpnts) / 2.])
        probsplt = array([0. * sum(gdat.maxmnumbpnts) / 2.])
        probmerg = array([0. * sum(gdat.maxmnumbpnts) / 2.])
        
        problgal = array([sum(gdat.maxmnumbpnts) / 2.])
        probbgal = array([sum(gdat.maxmnumbpnts) / 2.])
        probspec = array([sum(gdat.maxmnumbpnts) / 2.])
        probsind = array([sum(gdat.maxmnumbpnts) / 2.])
           
        gdat.probprop = concatenate((probmeanpnts, probfluxdistslop, probfluxdistbrek, probfluxdistsloplowr, probfluxdistslopuppr, \
                                                        probpsfipara, probnormback, probbrth, probdeth, probsplt, probmerg, problgal, probbgal, probspec, probsind))
        
        gdat.probprop /= sum(gdat.probprop)
       

def retr_randunitpsfipara(gdat):

    while True:
        randunitpsfipara = rand(gdat.numbpsfipara)
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
                    indx = m * gdat.numbpsfiparaevtt + i * gdat.numbformpara
                    thisbool = thisbool and randunitpsfipara[indx+indxpar1] > randunitpsfipara[indx+indxpar0]
            if thisbool:
                break

    return randunitpsfipara


def setp(gdat):
   
    # number of processes
    gdat.strgproc = os.uname()[1]
    if gdat.numbproc == None:
        if gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu':
            gdat.numbproc = 20
        else:
            gdat.numbproc = 1
    gdat.indxproc = arange(gdat.numbproc) 

    if gdat.exprtype == 'ferm':
        gdat.nameexpr = 'Fermi-LAT'
    if gdat.exprtype == 'sdss':
        gdat.nameexpr = 'SDSS'
    
    gdat.numbchrototl = 4
    gdat.numbchrollik = 7

    # half size of the spatial prior
    gdat.maxmgangfudi = gdat.maxmgang * (2. - gdat.margfact)
    gdat.maxmgangmarg = gdat.maxmgang * gdat.margfact

    # temp
    gdat.boolintpanglcosi = False

    # number of bins
    gdat.numbbins = 10

    gdat.minmnumbpnts = 1

    # the normalized offset for text annotation of point sources in the frames
    gdat.offstext = gdat.maxmgang * 0.05
    
    gdat.indxback = arange(gdat.numbback)
    
    gdat.numbevtt = gdat.indxevttincl.size
    gdat.indxevtt = arange(gdat.numbevtt)
    
    # axes
    ## longitude
    gdat.numblgal = 10
    gdat.binslgal = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numblgal + 1)

    ## latitude
    gdat.numbbgal = 10
    gdat.binsbgal = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbbgal + 1)

    ## radial
    gdat.numbgang = 10
    gdat.binsgang = linspace(-gdat.maxmgang, gdat.maxmgang, gdat.numbgang + 1)

    ## azimuthal
    gdat.numbaang = 10
    gdat.binsaang = linspace(0., 2. * pi, gdat.numbaang + 1)

    ## flux
    gdat.numbflux = 10
    gdat.indxflux = arange(gdat.numbflux)
    # temp
    gdat.fluxbinstype = 'logt'
    if gdat.fluxbinstype == 'logt':
        gdat.binsflux = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbflux + 1)
    if gdat.fluxbinstype == 'log2':
        gdat.binsflux = 1e-100 * 10.**(10.**(linspace(log10(log10(1e100 * gdat.minmflux)), log10(log10(1e100 * gdat.maxmflux)), gdat.numbflux + 1)))

    gdat.meanflux = sqrt(gdat.binsflux[1:] * gdat.binsflux[:-1])
    gdat.diffflux = gdat.binsflux[1:] - gdat.binsflux[:-1]
    ### pivot flux bin
    gdat.indxfluxpivt = gdat.numbflux / 2
    gdat.pivtflux = gdat.meanflux[gdat.indxfluxpivt]

    ## color
    gdat.numbsind = 10
    gdat.binssind = linspace(gdat.minmsind, gdat.maxmsind, gdat.numbsind + 1)
    gdat.meansind = (gdat.binssind[1:] + gdat.binssind[:-1]) / 2.
    gdat.diffsind = gdat.binssind[1:] - gdat.binssind[:-1]

    ## energy
    gdat.numbener = gdat.indxenerincl.size
    gdat.indxenerinclbins = empty(gdat.numbener+1, dtype=int)
    gdat.indxenerinclbins[0:-1] = gdat.indxenerincl
    gdat.indxenerinclbins[-1] = gdat.indxenerincl[-1] + 1
    gdat.binsener = gdat.binsenerfull[gdat.indxenerinclbins]
    gdat.diffener = (roll(gdat.binsener, -1) - gdat.binsener)[0:-1]

    gdat.meanener = sqrt(roll(gdat.binsener, -1) * gdat.binsener)[0:-1]
    gdat.indxener = arange(gdat.numbener, dtype=int)
    
    gdat.indxenerfluxdist = array([gdat.numbener / 2])
    gdat.enerfluxdist = gdat.meanener[gdat.indxenerfluxdist]
        
    factener = (gdat.meanener[gdat.indxenerfluxdist] / gdat.meanener)**2

    gdat.minmspec = gdat.minmflux * factener
    gdat.maxmspec = gdat.maxmflux * factener
    gdat.binsspec = gdat.binsflux[None, :] * factener[:, None]
    gdat.meanspec = empty((gdat.numbener, gdat.numbflux))
    for i in gdat.indxener:
        gdat.meanspec[i, :] = sqrt(gdat.binsspec[i, 1:] * gdat.binsspec[i, :-1])

    # temp
    if gdat.exprtype == 'sdss':
        gdat.diffener = ones(gdat.numbener)
    
    # angular gdat.deviation
    # temp -- check that gdat.numbangl does not degrade the performance
    gdat.numbangl = 1000
    gdat.binsangl = linspace(0., gdat.maxmangl, gdat.numbangl) # [rad]
    gdat.binsanglcosi = sort(cos(gdat.binsangl))
    
    # spatial priors
    gdat.minmlgal = -gdat.maxmgang
    gdat.maxmlgal = gdat.maxmgang
    gdat.minmbgal = -gdat.maxmgang
    gdat.maxmbgal = gdat.maxmgang
   
    # input data
    if gdat.datatype == 'inpt':
        
        path = gdat.pathdata + '/' + gdat.strgexpr
        gdat.exprflux = pf.getdata(path)
        
        if gdat.pixltype == 'heal':
            if gdat.exprflux.ndim != 3:
                print 'exprflux should be a 3D numpy array if pixelization is HealPix.'
                return
        else:
            if gdat.exprflux.ndim != 4:
                print 'exprflux should be a 4D numpy array if pixelization is Cartesian.'
                return
        
        if gdat.exprflux.ndim == 3:
            gdat.numbpixlheal = gdat.exprflux.shape[1]
            gdat.numbsideheal = int(sqrt(gdat.numbpixlheal / 12))
        elif gdat.exprflux.ndim == 4:
            gdat.numbsidecart = gdat.exprflux.shape[1]
            gdat.exprflux = gdat.exprflux.reshape((gdat.exprflux.shape[0], gdat.numbsidecart**2, gdat.exprflux.shape[3]))
            gdat.numbpixlheal = None
        else:
            print 'Input count map needs to be three (HealPix) or four (Cartesian) dimensional.'
            print 'Quitting...'
            return

        gdat.indxenerinclfull = arange(gdat.exprflux.shape[0])
        gdat.indxevttinclfull = arange(gdat.exprflux.shape[2])
        
        if gdat.pixltype == 'heal':
            gdat.numbsideheal = gdat.numbsideheal
        if gdat.pixltype == 'cart':
            gdat.numbsidecart = gdat.numbsidecart
        
    if gdat.datatype == 'mock':
        if gdat.exprtype == 'ferm':
            gdat.indxenerinclfull = arange(5)
            gdat.indxevttinclfull = arange(4)
        if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
            gdat.indxevttinclfull = arange(1)
            if dat.exprtype == 'sdss':
                gdat.indxenerinclfull = arange(3)
            else:
                gdat.indxenerinclfull = arange(2)

    # pizelization
    if gdat.pixltype == 'heal':
        gdat.numbpixlheal = gdat.numbsideheal**2 * 12
        gdat.apix = 4. * pi / gdat.numbpixlheal
    if gdat.pixltype == 'cart':
        gdat.binslgalcart = linspace(gdat.minmlgal, gdat.maxmlgal, gdat.numbsidecart + 1)
        gdat.binsbgalcart = linspace(gdat.minmbgal, gdat.maxmbgal, gdat.numbsidecart + 1)
        gdat.lgalcart = (gdat.binslgalcart[0:-1] + gdat.binslgalcart[1:]) / 2.
        gdat.bgalcart = (gdat.binsbgalcart[0:-1] + gdat.binsbgalcart[1:]) / 2.
        gdat.apix = (2. * gdat.maxmgang / gdat.numbsidecart)**2
        
    if gdat.pixltype == 'heal':
        
        lgalheal, bgalheal, gdat.numbpixlheal, gdat.apix = tdpy.util.retr_healgrid(gdat.numbsideheal)
        lgalheal = deg2rad(lgalheal)
        bgalheal = deg2rad(bgalheal)
   
        gdat.indxpixlrofi = where((fabs(lgalheal) < gdat.maxmgang) & (fabs(bgalheal) < gdat.maxmgang))[0]
        
        gdat.indxpixlrofimarg = where((fabs(lgalheal) < gdat.maxmgangmarg) & (fabs(bgalheal) < gdat.maxmgangmarg))[0]

        gdat.lgalgrid = lgalheal
        gdat.bgalgrid = bgalheal

    else:
        gdat.indxpixlrofi = arange(gdat.numbsidecart**2)
        indxsidecart = arange(gdat.numbsidecart)
        temp = meshgrid(indxsidecart, indxsidecart, indexing='ij')
        gdat.bgalgrid = gdat.bgalcart[temp[1].flatten()]
        gdat.lgalgrid = gdat.lgalcart[temp[0].flatten()]
    
    # plotting
    ## figure size
    gdat.plotsize = 7
    ## text
    if gdat.datatype == 'mock':
        gdat.truelabl = 'Mock'
    if gdat.datatype == 'inpt':
        if gdat.exprtype == 'ferm':
            gdat.truelabl = '3FGL'
        if gdat.exprtype == 'sdss':
            gdat.truelabl = 'Hubble'
        if gdat.exprtype == 'chan':
            gdat.truelabl = 'Chandra'

    if gdat.exprtype == 'ferm':
        gdat.strganglunit = '$^o$'
    if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
        gdat.strganglunit = '$^{\prime\prime}$'
            
    if gdat.exprtype == 'ferm':
        gdat.enerfact = 1.
    if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
        gdat.enerfact = 1e3
    
    if gdat.strganglunit == '$^o$':
        gdat.anglfact = 180. / pi
    if gdat.strganglunit == '$^{\prime\prime}$':
        gdat.anglfact = 3600 * 180. / pi

    gdat.binsanglplot = gdat.anglfact * gdat.binsangl

    if gdat.regitype == 'igal':
        gdat.longlabl = '$l$'
        gdat.latilabl = '$b$'
    if gdat.regitype == 'ngal':
        gdat.longlabl = r'$\nu$'
        gdat.latilabl = r'$\mu$'
    
    gdat.longlabl += ' [%s]' % gdat.strganglunit
    gdat.latilabl += ' [%s]' % gdat.strganglunit
    
    gdat.truelablvari = gdat.truelabl + ' variable'
    gdat.truelablmiss = gdat.truelabl + ' miss'
    gdat.truelablbias = gdat.truelabl + ' off'
    gdat.truelablhits = gdat.truelabl + ' hit'
    gdat.truelablmult = gdat.truelabl + ' mult'

    # minimum angular distance from the center of the ROI
    gdat.minmgang = 1e-3

    # energy bin string
    gdat.enerstrg, gdat.binsenerstrg = retr_enerstrg(gdat)
    
    # PSF class string
    gdat.evttstrg = []
    for m in gdat.indxevtt:
        gdat.evttstrg.append('PSF%d' % gdat.indxevttincl[m])
        
    if gdat.exprtype == 'ferm':
        retr_fermpsfn(gdat)
    if gdat.exprtype == 'sdss':
        retr_sdsspsfn(gdat)

    gdat.spltjcbnfact = log(2.**(2 - gdat.numbener))
    
    # number of components
    # temp -- 3->4 4->5
    gdat.numbcomp = 4 + gdat.numbener
    gdat.numbcompcolr = 5
    gdat.jcbnsplt = 2.**(2 - gdat.numbener)
    
    # component indices
    gdat.indxcomplgal = 0
    gdat.indxcompbgal = 1
    gdat.indxcompspec = 2 + gdat.indxener
    gdat.indxcompflux = 2 + gdat.indxenerfluxdist
    gdat.indxcompsind = 2 + gdat.numbener
    
    # population index vector
    gdat.indxpopl = arange(gdat.numbpopl, dtype=int)
    if gdat.datatype == 'mock':
        gdat.mockindxpopl = arange(gdat.mocknumbpopl, dtype=int)

    # convenience factors for CDF and ICDF transforms
    ## mean number of PS
    gdat.factmeanpnts = log(gdat.maxmmeanpnts / gdat.minmmeanpnts)
    
    ## background normalization
    gdat.factnormback = log(gdat.maxmnormback / gdat.minmnormback)
    
    ## PS parameters
    gdat.factgang = log(gdat.maxmgangmarg / gdat.minmgang)

    gdat.factfluxdistslop = arctan(gdat.maxmfluxdistslop) - arctan(gdat.minmfluxdistslop)
    gdat.factfluxdistbrek = log(gdat.maxmfluxdistbrek / gdat.minmfluxdistbrek)
    gdat.factfluxdistsloplowr = arctan(gdat.maxmfluxdistsloplowr) - arctan(gdat.minmfluxdistsloplowr)
    gdat.factfluxdistslopuppr = arctan(gdat.maxmfluxdistslopuppr) - arctan(gdat.minmfluxdistslopuppr)
    gdat.factsind = arctan(gdat.maxmsind) - arctan(gdat.minmsind)
    
    if gdat.datatype == 'mock':
        gdat.mocksindcdfnnormminm = 0.5 * (sp.special.erf((gdat.minmsind - gdat.mocksinddistmean) / gdat.mocksinddiststdv / sqrt(2.)) + 1.)
        gdat.mocksindcdfnnormmaxm = 0.5 * (sp.special.erf((gdat.maxmsind - gdat.mocksinddistmean) / gdat.mocksinddiststdv / sqrt(2.)) + 1.)
        gdat.mocksindcdfnnormdiff = gdat.mocksindcdfnnormmaxm - gdat.mocksindcdfnnormminm
    
    # temp -- these must be updated when updating hyperparameters on the color distribution
    gdat.sindcdfnnormminm = 0.5 * (sp.special.erf((gdat.minmsind - gdat.sinddistmean) / gdat.sinddiststdv / sqrt(2.)) + 1.)
    gdat.sindcdfnnormmaxm = 0.5 * (sp.special.erf((gdat.maxmsind - gdat.sinddistmean) / gdat.sinddiststdv / sqrt(2.)) + 1.)
    gdat.sindcdfnnormdiff = gdat.sindcdfnnormmaxm - gdat.sindcdfnnormminm

    # construct the mock PSF
    if gdat.datatype == 'mock':

        # mock PSF parameters
        if gdat.mockpsfntype == 'singgaus':
            gdat.numbmockformpara = 1
        elif gdat.mockpsfntype == 'singking':
            gdat.numbmockformpara = 2 
        elif gdat.mockpsfntype == 'doubgaus':
            gdat.numbmockformpara = 3
        elif gdat.mockpsfntype == 'gausking':
            gdat.numbmockformpara = 4
        elif gdat.mockpsfntype == 'doubking':
            gdat.numbmockformpara = 5

        gdat.numbmockpsfiparaevtt = gdat.numbener * gdat.numbmockformpara
        gdat.numbmockpsfipara = gdat.numbmockpsfiparaevtt * gdat.numbevtt
        gdat.indxmockpsfipara = arange(gdat.numbmockpsfipara)   

        if gdat.exprtype == 'ferm':
            gdat.mockpsfipara = gdat.fermpsfipara
        if gdat.exprtype == 'sdss':
            gdat.mockpsfipara = gdat.sdsspsfipara
    
        if gdat.pixltype == 'heal':
            gdat.numbpixlfull = 12 * gdat.numbsideheal**2
        if gdat.pixltype == 'cart':
            gdat.numbpixlfull = gdat.numbsidecart**2
    if gdat.datatype == 'inpt':
        gdat.numbpixlfull = gdat.exprflux.shape[1]

    # exposure
    if gdat.strgexpo == 'unit':
        if gdat.datatype == 'mock':
            if gdat.pixltype == 'heal':
                gdat.expo= ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
            if gdat.pixltype == 'cart':
                gdat.expo = ones((gdat.numbenerfull, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevttfull))
        if gdat.datatype == 'inpt':
            gdat.expo = ones_like(gdat.exprflux)
    else:
        path = gdat.pathdata + '/' + gdat.strgexpo
        gdat.expo = pf.getdata(path)
        if amin(gdat.expo) == amax(gdat.expo):
            print 'Bad input exposure map.'
            return
        if gdat.pixltype == 'cart':
            gdat.expo = gdat.expo.reshape((gdat.expo.shape[0], -1, gdat.expo.shape[-1]))

    # backgrounds
    gdat.backflux = []
    for c in gdat.indxback:
        if gdat.strgback[c] == 'unit':
            if gdat.datatype == 'mock':
                if gdat.datatype == 'heal':
                    backfluxtemp = ones((gdat.numbenerfull, gdat.numbpixlfull, gdat.numbevttfull))
                if gdat.datatype == 'cart':
                    backfluxtemp = ones((gdat.numbenerfull, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevttfull))
            if gdat.datatype == 'inpt':
                backfluxtemp = ones_like(gdat.exprflux)
        else:
            path = gdat.pathdata + '/' + gdat.strgback[c]
            backfluxtemp = pf.getdata(path)
            if gdat.pixltype == 'cart':
                backfluxtemp = backfluxtemp.reshape((backfluxtemp.shape[0], -1, backfluxtemp.shape[-1]))
        gdat.backflux.append(backfluxtemp)
    
    # only include desired energy and PSF class bins 
    gdat.indxcubeincl = meshgrid(gdat.indxenerincl, arange(gdat.expo.shape[1]), gdat.indxevttincl, indexing='ij')
    if gdat.datatype == 'inpt':
        gdat.exprflux = gdat.exprflux[gdat.indxcubeincl]
    gdat.expo = gdat.expo[gdat.indxcubeincl]
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcubeincl]
    
    # exclude voxels with vanishing exposure
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
        path = os.environ["PCAT_DATA_PATH"] + '/pixlcnvt_%09g.p' % (gdat.maxmgang)

        if os.path.isfile(path):
            fobj = open(path, 'rb')
            gdat.pixlcnvt = cPickle.load(fobj)
            fobj.close()
        else:
            gdat.pixlcnvt = zeros(gdat.numbpixlheal, dtype=int) - 1
            numbpixlmarg = gdat.indxpixlrofimarg.size

            for k in range(numbpixlmarg):
                dist = retr_angldistunit(gdat, lgalheal[gdat.indxpixlrofimarg[k]], bgalheal[gdat.indxpixlrofimarg[k]], gdat.indxpixl)
                gdat.pixlcnvt[gdat.indxpixlrofimarg[k]] = argmin(dist)
            fobj = open(path, 'wb')
            cPickle.dump(gdat.pixlcnvt, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
            fobj.close()
   
    
    print 'hey'
    print retr_indxpixl(gdat, 0., 0.)
    print retr_indxpixl(gdat, deg2rad(10.), deg2rad(10.))
    print retr_indxpixl(gdat, deg2rad(20.), deg2rad(20.))
    print retr_indxpixl(gdat, deg2rad(30.), deg2rad(30.))

    if gdat.datatype == 'inpt':
        # temp
        #gdat.datacnts = gdat.datacnts[gdat.indxcuberofi]
        gdat.exprflux = gdat.exprflux[gdat.indxcuberofi]
    
    gdat.expofull = copy(gdat.expo)
    gdat.expo = gdat.expo[gdat.indxcuberofi]
    for c in gdat.indxback:
        gdat.backflux[c] = gdat.backflux[c][gdat.indxcuberofi]

    # get the experimental catalog
    if gdat.exprinfo:
        if gdat.exprtype == 'ferm':
            retr_fermdata(gdat)
        gdat.exprgang = retr_gang(gdat.exprlgal, gdat.exprbgal)
        gdat.expraang = retr_aang(gdat.exprlgal, gdat.exprbgal)

    # flag to indicate whether information from a deterministic catalog will be used or not
    gdat.trueinfo = gdat.datatype == 'mock' or gdat.exprinfo
    
    # get count data
    ## input data
    if gdat.datatype == 'inpt':
        gdat.datacnts = gdat.exprflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]
    
    ## mock data
    if gdat.datatype == 'mock':

        if gdat.mocknumbpnts == None:
            gdat.mocknumbpnts = empty(gdat.numbpopl)
            for l in gdat.indxpopl:
                gdat.mocknumbpnts[l] = random_integers(gdat.minmnumbpnts, gdat.maxmnumbpnts[l])
            
        gdat.truemeanpnts = gdat.mocknumbpnts
        
        # if mock FDF is not specified by the user, randomly seed it from the prior
        # temp -- make this section compatible with mock variables being None
        for l in gdat.indxpopl:
            if gdat.mockfluxdisttype[l] == 'powr':
                if gdat.mockfluxdistslop[l] == None:
                    gdat.mockfluxdistslop[l] = icdf_atan(rand(), gdat.minmfluxdistslop[l], gdat.factfluxdistslop[l])
            if gdat.mockfluxdisttype[l] == 'brok':
                if gdat.mockfluxdistbrek[l] == None:
                    gdat.mockfluxdistbrek[l] = icdf_atan(rand(), gdat.minmfluxdistbrek[l], gdat.factfluxdistbrek[l])
                if gdat.mockfluxdistsloplowr[l] == None:
                    gdat.mockfluxdistsloplowr[l] = icdf_atan(rand(), gdat.minmfluxdistsloplowr[l], gdat.factfluxdistsloplowr[l])
                if gdat.mockfluxdistslopuppr[l] == None:
                    gdat.mockfluxdistslopuppr[l] = icdf_atan(rand(), gdat.minmfluxdistslopuppr[l], gdat.factfluxdistslopuppr[l])

        gdat.truefluxdistslop = gdat.mockfluxdistslop
        gdat.truefluxdistbrek = gdat.mockfluxdistbrek
        gdat.truefluxdistsloplowr = gdat.mockfluxdistsloplowr
        gdat.truefluxdistslopuppr = gdat.mockfluxdistslopuppr
    
        if gdat.mockpsfipara == None: 
            gdat.mockpsfntype = psfntpye
            numbmockpsfipara = gdat.numbpsfipara
            gdat.mockpsfipara = empty(numbmockpsfipara)
            for k in arange(numbmockpsfipara):
                gdat.mockpsfipara[k] = icdf_psfipara(gdat, rand(), k)
   
        if gdat.mocknormback == None:
            for c in gdat.indxback:
                gdat.mocknormback[c, :] = icdf_logt(rand(gdat.numbener), gdat.minmnormback[c], gdat.factnormback[c])

        mockcnts = [[] for l in gdat.indxpopl]
        mocklgal = [[] for l in gdat.indxpopl]
        mockbgal = [[] for l in gdat.indxpopl]
        mockgang = [[] for l in gdat.indxpopl]
        mockaang = [[] for l in gdat.indxpopl]
        mockspec = [[] for l in gdat.indxpopl]
        mocksind = [[] for l in gdat.indxpopl]
        for l in gdat.mockindxpopl:
            if gdat.mockspatdisttype[l] == 'unif':
                mocklgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg)
                mockbgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
            if gdat.mockspatdisttype[l] == 'disc':
                mockbgal[l] = icdf_logt(rand(gdat.mocknumbpnts[l]), gdat.minmgang, gdat.factgang) * choice(array([1., -1.]), size=gdat.mocknumbpnts[l])
                mocklgal[l] = icdf_self(rand(gdat.mocknumbpnts[l]), -gdat.maxmgangmarg, 2. * gdat.maxmgangmarg) 
            if gdat.mockspatdisttype[l] == 'gang':
                mockgang[l] = icdf_logt(rand(gdat.mocknumbpnts[l]), gdat.minmgang, gdat.factgang)
                mockaang[l] = icdf_self(rand(gdat.mocknumbpnts[l]), 0., 2. * pi)
                mocklgal[l], mockbgal[l] = retr_lgalbgal(mockgang[l], mockaang[l])

            mockspec[l] = empty((gdat.numbener, gdat.mocknumbpnts[l]))
            if gdat.mockfluxdisttype[l] == 'powr':
                mockspec[l][gdat.indxenerfluxdist[0], :] = icdf_flux_powr(gdat, rand(gdat.mocknumbpnts[l]), gdat.mockfluxdistslop[l])
            if gdat.mockfluxdisttype[l] == 'brok':
                mockspec[l][gdat.indxenerfluxdist[0], :] = icdf_flux_brok(gdat, rand(gdat.mocknumbpnts[l]), gdat.mockfluxdistbrek[l], \
                                                                                                gdat.mockfluxdistsloplowr[l], gdat.mockfluxdistslopuppr[l])
            mocksind[l] = icdf_eerr(rand(gdat.mocknumbpnts[l]), gdat.mocksinddistmean[l], gdat.mocksinddiststdv[l], gdat.mocksindcdfnnormminm[l], gdat.mocksindcdfnnormdiff[l])
        
            if gdat.verbtype > 1:
                print 'mocksind[l]'
                print mocksind[l]
                print 'mockspec[l]'
                print mockspec[l]
                print

            if gdat.mockspectype[l] == 'powr':
                mockspec[l] = retr_spec(gdat, mockspec[l][gdat.indxenerfluxdist[0], :], mocksind[l])
            if gdat.mockspectype[l] == 'curv':
                mockspec[l] = retr_speccurv(gdat, mockspec[l][gdat.indxenerfluxdist[0], :], mocksind[l], mockcurv[l])
            if gdat.mockspectype[l] == 'expo':
                mockspec[l] = retr_spec(gdat, mockspec[l][gdat.indxenerfluxdist[0], :], mocksind[l], mockbrek)[l]
            
            indxpixltemp = retr_indxpixl(gdat, mockbgal[l], mocklgal[l])
            
            if False:
                print 'hey'
                print 'mock'
                print 'mockbgal[l]'
                print mockbgal[l]
                print 'mocklgal[l]'
                print mocklgal[l]
                print 'indxpixltemp'
                print indxpixltemp
                print 

            mockcnts[l] = mockspec[l][:, :, None] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None]
        
        mockpntsflux = retr_pntsflux(gdat, concatenate(mocklgal), concatenate(mockbgal), concatenate(mockspec, axis=1), gdat.mockpsfipara, gdat.mockpsfntype)
        mocktotlflux = retr_rofi_flux(gdat, gdat.mocknormback, mockpntsflux, gdat.indxcube)
        mocktotlcnts = mocktotlflux * gdat.expo * gdat.apix * gdat.diffener[:, None, None] # [1]

        gdat.datacnts = zeros((gdat.numbener, gdat.numbpixl, gdat.numbevtt))
        for i in gdat.indxener:
            for k in range(gdat.numbpixl):
                for m in gdat.indxevtt:
                    gdat.datacnts[i, k, m] = poisson(mocktotlcnts[i, k, m])
                    
        if gdat.trueinfo:
            
            gdat.truelgal = []
            gdat.truebgal = []
            gdat.truesind = []
            for l in gdat.indxpopl:
                gdat.truelgal.append(mocklgal[l])
                gdat.truebgal.append(mockbgal[l])
                gdat.truesind.append(mocksind[l])
                    
            gdat.truestrg = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truestrgclss = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.truestrgassc = [array([None for n in range(gdat.mocknumbpnts[l])], dtype=object) for l in gdat.indxpopl]
            gdat.indxtruepntstimevari = [array([])] * gdat.numbpopl
            gdat.truenumbpnts = gdat.mocknumbpnts
            
            if gdat.mockfluxdisttype[l] == 'powr':
                gdat.truefluxdistslop = gdat.mockfluxdistslop
            if gdat.mockfluxdisttype[l] == 'brok':
                gdat.truefluxdistbrek = gdat.mockfluxdistbrek
                gdat.truefluxdistsloplowr = gdat.mockfluxdistsloplowr
                gdat.truefluxdistslopuppr = gdat.mockfluxdistslopuppr
            gdat.truecnts = mockcnts
               
            gdat.truespec = []
            for l in gdat.indxpopl:
                gdat.truespectemp = empty((3, gdat.numbener, gdat.mocknumbpnts[l]))
                gdat.truespectemp[:] = mockspec[l][None, :, :]
                gdat.truespec.append(gdat.truespectemp)
            
            gdat.truepsfipara = gdat.mockpsfipara
            gdat.truepsfntype = gdat.mockpsfntype
            gdat.truenormback = gdat.mocknormback
        
    if gdat.trueinfo:
        gdat.indxtruepntsfudi = []
        for l in gdat.indxpopl:
            indxtruepntstemp = where((fabs(gdat.truelgal[l]) < gdat.maxmgangfudi) & (fabs(gdat.truebgal[l]) < gdat.maxmgangfudi))[0]
            gdat.indxtruepntsfudi.append(indxtruepntstemp)

    # spatially averaged background flux 
    gdat.backfluxmean = zeros((gdat.numbback, gdat.numbener))
    for c in gdat.indxback:
        for i in gdat.indxener:
            gdat.backfluxmean[c, i] += sum(gdat.backflux[c][i, :, :] * gdat.expo[i, :, :]) / sum(gdat.expo[i, :, :])

    # maximum number of point sources that can be modified at once
    gdat.numbmodipnts = int(max(3, sum(gdat.maxmnumbpnts)))
    
    # construct the PSF model
    retr_psfimodl(gdat)

    # proposals
    retr_propmodl(gdat)
    
    # factors in the prior expression
    if gdat.priotype == 'expo':
        gdat.priofactexpo = gdat.numbcompcolr / 2.
    gdat.priofactmeanpnts = log(1. / (log(gdat.maxmmeanpnts) - log(gdat.minmmeanpnts)))
    gdat.priofactlgalbgal = 2. * log(1. / 2. / gdat.maxmgang)
    gdat.priofactfluxdistslop = gdat.numbener * log(1. / (arctan(gdat.maxmfluxdistslop) - arctan(gdat.minmfluxdistslop)))
    # temp -- brok terms are missing

    # initialize the counter
    cntr = tdpy.util.cntr()
    
    # sample vector indices  
    gdat.indxsampnumbpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampmeanpnts = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    
    # dummy definitions
    gdat.indxsampfluxdistslop = -1.
    gdat.indxsampfluxdistbrek = -1
    gdat.indxsampfluxdistsloplowr = -1
    gdat.indxsampfluxdistslopuppr = -1

    gdat.indxsampfluxdistslop = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistbrek = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistsloplowr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsampfluxdistslopuppr = arange(gdat.numbpopl) + cntr.incr(gdat.numbpopl)
    gdat.indxsamppsfipara = arange(gdat.numbpsfipara) + cntr.incr(gdat.numbpsfipara)
    gdat.indxsampnormback = arange(gdat.numbback * gdat.numbener).reshape((gdat.numbback, gdat.numbener)) + cntr.incr(gdat.numbback * gdat.numbener)

    gdat.maxmnumbcomp = gdat.maxmnumbpnts * gdat.numbcomp
    gdat.indxsampcompinit = amax(gdat.indxsampnormback) + 1
    
    # maximum number of parameters
    gdat.numbpara = int(gdat.indxsampcompinit + sum(gdat.maxmnumbcomp))
    gdat.indxpara = arange(gdat.numbpara)

    if gdat.numbburn == None:
        gdat.numbburn = min(200000, gdat.numbswep - 1)

    if gdat.factthin == None:
        gdat.factthin = min(2 * gdat.numbpara, gdat.numbswep - gdat.numbburn)

    # run tag
    gdat.rtag = retr_rtag(gdat, None)
    
    # plots
    if gdat.makeplot:
        if (gdat.strgproc == 'fink1.rc.fas.harvard.edu' or gdat.strgproc == 'fink2.rc.fas.harvard.edu') and getpass.getuser() == 'tansu':
            pathplotbase = '/n/pan/www/tansu/imag/pcat/'
        else:
            pathplotbase = gdat.pathdata + '/imag/'
        gdat.pathplot = pathplotbase + gdat.strgtime + '_' + gdat.strgcnfg + '_' + gdat.rtag + '/'
        cmnd = 'mkdir -p ' + gdat.pathplot
        os.system(cmnd)

    # number of samples to be saved
    gdat.numbsamp = (gdat.numbswep - gdat.numbburn) / gdat.factthin
    gdat.indxsamp = arange(gdat.numbsamp)
    gdat.numbsamptotl = gdat.numbsamp * gdat.numbproc
    gdat.indxsamptotl = arange(gdat.numbsamptotl)
    gdat.numbsweptotl = gdat.numbswep * gdat.numbproc

    # rescale the positional update scale
    gdat.stdvlbhl /= 2. * gdat.maxmgang

    # determine proposal probabilities
    gdat.probpropminm = copy(gdat.probprop)
    gdat.probpropmaxm = copy(gdat.probprop)
   
    ## do not allow death or merge when the number of PS is at its minimum or 
    ## births or splits when the number of PS is at its maximum
    gdat.probpropminm[[gdat.indxpropdeth, gdat.indxpropmerg]] = 0.
    gdat.probpropmaxm[[gdat.indxpropbrth, gdat.indxpropsplt]] = 0.
    gdat.probprop /= sum(gdat.probprop)
    gdat.probpropmaxm /= sum(gdat.probpropmaxm)
    gdat.probpropminm /= sum(gdat.probpropminm)
    
    gdat.minmlgalmarg = -gdat.maxmgangmarg
    gdat.maxmlgalmarg = gdat.maxmgangmarg
    gdat.minmbgalmarg = -gdat.maxmgangmarg
    gdat.maxmbgalmarg = gdat.maxmgangmarg

    # temp
    gdat.tracsamp = False
    
    # center of the ROI
    if gdat.regitype == 'igal':
        gdat.cntrlghp, gdat.cntrbghp = 0., 0.
    else:
        gdat.cntrlghp, gdat.cntrbghp = 0., 90.
   
    # plot settings
    ## marker opacity
    gdat.mrkralph = 0.5
    ## marker size
    gdat.minmmrkrsize = 50
    gdat.maxmmrkrsize = 250
    ## ROI
    gdat.exttrofi = array([gdat.minmlgal, gdat.maxmlgal, gdat.minmbgal, gdat.maxmbgal])
    gdat.exttrofi *= gdat.anglfact 
    gdat.frambndr = gdat.maxmgang * gdat.anglfact
    gdat.frambndrmarg = gdat.maxmgangmarg * gdat.anglfact
    
    # convenience variables
    gdat.fluxhistmodl = empty(gdat.numbflux)
    
    if gdat.exprtype == 'ferm':
        gdat.numbfluxprox = 3
    if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
        gdat.numbfluxprox = 1
    gdat.indxfluxprox = arange(gdat.numbfluxprox)
    gdat.binsfluxprox = logspace(log10(gdat.minmflux), log10(gdat.maxmflux), gdat.numbfluxprox + 1)
    gdat.meanfluxprox = sqrt(gdat.binsfluxprox[1:] * gdat.binsfluxprox[:-1])
    
    gdat.maxmangleval = empty(gdat.numbfluxprox)
    for h in gdat.indxfluxprox:
        if gdat.specfraceval == 0:
            gdat.maxmangleval[h] = 3. * gdat.maxmgangmarg
        else:
            if gdat.exprtype == 'ferm':
                frac = gdat.specfraceval * gdat.binsfluxprox[0] / gdat.binsfluxprox[h+1]
                psfnwdth = retr_psfnwdth(gdat, gdat.fermpsfn, frac)
                gdat.indxmaxmangl = unravel_index(argmax(psfnwdth), psfnwdth.shape)
                gdat.maxmangleval[h] = psfnwdth[gdat.indxmaxmangl]
            if gdat.exprtype == 'sdss' or gdat.exprtype == 'chan':
                gdat.maxmangleval[h] = 10. / gdat.anglfact

    # pixels whose posterior predicted emission will be saved
    gdat.numbpixlsave = min(1000, gdat.numbpixl)
    gdat.indxpixlsave = choice(arange(gdat.numbpixlsave), size=gdat.numbpixlsave)

    factener = (gdat.meanener / gdat.meanener[gdat.indxenerfluxdist[0]])**(-2.)
    
    gdat.minmcnts = 1e-1 * factener
    gdat.maxmcnts = 1e4 * factener
    # temp 
    gdat.minmcnts = gdat.minmflux * 1e-2 * mean(mean(gdat.expo, 1), 1) * gdat.diffener * factener
    gdat.maxmcnts = gdat.maxmflux * 1e2 * mean(mean(gdat.expo, 1), 1) * gdat.diffener * factener

    gdat.binscnts = zeros((gdat.numbener, gdat.numbflux + 1))
    for i in gdat.indxener:
        gdat.binscnts[i, :] = logspace(log10(gdat.minmcnts[i]), log10(gdat.maxmcnts[i]), gdat.numbflux + 1) # [1]
       
    ## Real data
    # true data
    if gdat.trueinfo:
        if gdat.datatype == 'inpt':
            gdat.truenumbpnts = None
            gdat.truemeanpnts = None
            gdat.truenormback = None
    
            if gdat.exprtype == 'ferm':
                gdat.truenumbpnts = array([gdat.exprnumbpnts], dtype=int)
                gdat.truelgal = [gdat.exprlgal]
                gdat.truebgal = [gdat.exprbgal]
                gdat.truespec = [gdat.exprspec]
                gdat.truecnts = [gdat.exprcnts]
                gdat.truesind = [gdat.exprsind]
                gdat.truestrg = [gdat.exprstrg]
                gdat.truestrgclss = [gdat.exprstrgclss]
                gdat.truestrgassc = [gdat.exprstrgassc]
                gdat.indxtruevari = [gdat.indxexprvari]
                gdat.truepsfipara = gdat.fermpsfipara
                gdat.truepsfntype = 'doubking'
                
        if gdat.datatype == 'mock':
            gdat.truepsfn = retr_psfn(gdat, gdat.truepsfipara, gdat.indxener, gdat.binsangl, gdat.mockpsfntype)
        else:
            if gdat.exprtype == 'ferm':
                gdat.truepsfn = gdat.fermpsfn
            elif gdat.exprtype == 'sdss':
                gdat.truepsfn = gdat.sdsspsfn
            
        if gdat.truepsfn != None:
            gdat.truefwhm = 2. * retr_psfnwdth(gdat, gdat.truepsfn, 0.5)
        
        if gdat.trueinfo:
            truebackcnts = []
            gdat.truesigm = []
            for l in gdat.indxpopl:
                indxpixltemp = retr_indxpixl(gdat, gdat.truebgal[l], gdat.truelgal[l])
                truebackcntstemp = zeros((gdat.numbener, gdat.truenumbpnts[l], gdat.numbevtt))
                for c in gdat.indxback:
                    truebackcntstemp += gdat.backflux[c][:, indxpixltemp, :] * gdat.expo[:, indxpixltemp, :] * gdat.diffener[:, None, None] * pi * \
                                                                                                                                gdat.truefwhm[:, None, :]**2 / 4.
                truebackcnts.append(truebackcntstemp)
                gdat.truesigm.append(gdat.truecnts[l] / sqrt(truebackcntstemp))
    else:
        gdat.truepsfn = None
      
    # sanity checks
    # temp
    if (fabs(gdat.datacnts - rint(gdat.datacnts)) > 1e-3).any():
        print 'Fractional counts!'

    if amin(gdat.datacnts) < 0.:
        print 'Negative counts!'

    gdat.datafluxmean = sum(sum(gdat.datacnts, 1), 1) / sum(sum(gdat.expo, 1), 1) / gdat.apix / gdat.diffener
    gdat.datacntsmean = mean(sum(gdat.datacnts, 2), 1)
    gdat.datacntssatu = ceil((amax(sum(gdat.datacnts, 2), 1) - gdat.datacntsmean) * 0.05 + gdat.datacntsmean)
    gdat.resicntssatu = ceil(gdat.datacntssatu * 0.2)
    
    # auxiliary variables for plots
    # temp
    if False:
        if gdat.pixltype == 'heal':
            for i in gdat.indxener:
                for m in gdat.indxevtt:
                    gdat.datacntscarttemp = tdpy.util.retr_cart(gdat.datacnts[i, :, m], gdat.indxpixlrofi, numbsideinpt=gdat.numbsideheal, \
                        minmlgal=gdat.anglfact*gdat.minmlgal, maxmlgal=gdat.anglfact*gdat.maxmlgal, minmbgal=gdat.anglfact*gdat.minmbgal, maxmbgal=gdat.anglfact*gdat.maxmbgal)
                    if i == 0 and m == 0:
                        gdat.datacntscart = zeros((gdat.datacntscarttemp.shape[0], gdat.datacntscarttemp.shape[1], gdat.numbener, gdat.numbevtt))
                    gdat.datacntscart[:, :, i, m] = gdat.datacntscarttemp
        if gdat.pixltype == 'cart':
            gdat.datacntscart = zeros((gdat.numbener, gdat.numbpixlfull, gdat.numbevtt)) 
            gdat.datacntscart[gdat.indxcuberofi] = gdat.datacnts
            gdat.datacntscart = gdat.datacntscart.reshape((gdat.numbener, gdat.numbsidecart, gdat.numbsidecart, gdat.numbevtt))
            gdat.datacntscart = swapaxes(swapaxes(gdat.datacntscart, 0, 2), 0, 1)

        for i in gdat.indxener:
            indxdatacntscartsatu = where(gdat.datacntscart[:, :, i, :] > gdat.datacntssatu[i])
        gdat.datacntscart[indxdatacntscartsatu[0], indxdatacntscartsatu[1], i, indxdatacntscartsatu[2]] = gdat.datacntssatu[i]

    # make a look-up table of nearby pixels for each pixel
    if gdat.specfraceval == 0:
        path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%06d_%s.p' % (gdat.numbpixl, gdat.pixltype)
    else:
        path = os.environ["PCAT_DATA_PATH"] + '/indxpixlprox_%06d_%s_%.7g_%.7g_%02d.p' % (gdat.numbpixl, gdat.pixltype, gdat.minmflux, gdat.maxmflux, gdat.numbfluxprox)

    global indxpixlprox
    if os.path.isfile(path):
        print 'Retrieving previously computed pixel look-up table...'
        fobj = open(path, 'rb')
        indxpixlprox = cPickle.load(fobj)
        fobj.close()
    else:
        print 'Computing the look-up table...'
        indxpixlprox = [[] for h in range(gdat.numbfluxprox)]
        for j in gdat.indxpixl:
            dist = retr_angldistunit(gdat, gdat.lgalgrid[j], gdat.bgalgrid[j], gdat.indxpixl)
            dist[j] = 0.
            for h in range(gdat.numbfluxprox):
                indxpixlproxtemp = where(dist < gdat.maxmangleval[h])[0]
                indxpixlproxtemp = indxpixlproxtemp[argsort(dist[indxpixlproxtemp])]
                indxpixlprox[h].append(indxpixlproxtemp)
        fobj = open(path, 'wb')
        cPickle.dump(indxpixlprox, fobj, protocol=cPickle.HIGHEST_PROTOCOL)
        fobj.close()
        
    if gdat.verbtype > 1:
        print 'Memory budget: indxpixlprox'
        totl = 0.
        for h in gdat.indxfluxprox:
            for n in gdat.indxpixl:
                totl += sys.getsizeof(indxpixlprox[h][n]) / 2.**20
        print '%.4g MB' % totl


def init_fram(gdat, gdatmodi, indxevttplot, indxenerplot, strgplot, indxpoplplot=None):

    figr, axis = plt.subplots(figsize=(1.3 * gdat.plotsize, 1.3 * gdat.plotsize))
    axis.set_xlabel(gdat.longlabl)
    axis.set_ylabel(gdat.latilabl)
    axis.set_xlim([gdat.frambndrmarg, -gdat.frambndrmarg])
    axis.set_ylim([-gdat.frambndrmarg, gdat.frambndrmarg])
    if indxevttplot == None:
        if indxenerplot == None:
            axis.set_title('')
        else:
            axis.set_title(gdat.binsenerstrg[indxenerplot])
    else:
        titl = gdat.binsenerstrg[indxenerplot]
        if gdat.exprtype == 'ferm':
            titl += ', ' + gdat.evttstrg[indxevttplot]
        axis.set_title(titl)

    axis.axvline(gdat.frambndr, ls='--', alpha=gdat.mrkralph, color='black')
    axis.axvline(-gdat.frambndr, ls='--', alpha=gdat.mrkralph, color='black')
    axis.axhline(gdat.frambndr, ls='--', alpha=gdat.mrkralph, color='black')
    axis.axhline(-gdat.frambndr, ls='--', alpha=gdat.mrkralph, color='black')

    if indxpoplplot == None:
        strg = ''
    else:
        strg = 'pop%d' % indxpoplplot
    if indxevttplot == None:
        path = gdat.pathplot + strgplot + strg + '_%dA_' % (gdat.indxenerincl[indxenerplot]) + gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    else:
        path = gdat.pathplot + strgplot + strg + '_%d%d_' % (gdat.indxenerincl[indxenerplot], gdat.indxevttincl[indxevttplot]) + \
                                                                                                                    gdat.rtag + '_%09d.pdf' % gdatmodi.cntrswep
    
    return figr, axis, path


def supr_fram(gdat, gdatmodi, axis, indxenerplot, indxpoplplot):

    # temp
    if True:
        if gdat.trueinfo:
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label='Sample', marker='+', linewidth=2, color='b')
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
            axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablhits, marker='*', linewidth=2, color='g')
            # temp
            #if gdat.indxtruepntstimevari[indxpoplplot].size > 0:
            #    axis.scatter(gdat.maxmgang * 5., gdat.maxmgang * 5, s=50, alpha=gdat.mrkralph, label=gdat.truelablvari, marker='*', linewidth=2, color='y')
        axis.legend(bbox_to_anchor=[0.5, 1.1], loc='center', ncol=2)
        
    # true catalog
    if gdat.trueinfo:
        ## get the true catalog
        mrkrsize = retr_mrkrsize(gdat, gdat.truespec[indxpoplplot][0, gdat.indxenerfluxdist, :].flatten())
        lgal = copy(gdat.truelgal[indxpoplplot])
        bgal = copy(gdat.truebgal[indxpoplplot])
        numbpnts = int(gdat.truenumbpnts[indxpoplplot])
        
        ## associations
        ### missed
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].miss
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.mrkralph, label=gdat.truelablmiss, marker='x', linewidth=2, color='g')
        
        ### biased
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].bias[indxenerplot]
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.mrkralph, \
                                                                                    label=gdat.truelablbias, marker='o', linewidth=2, color='g', facecolor='none')
        
        ### hit
        indx = gdatmodi.indxtruepntsassc[indxpoplplot].hits[indxenerplot]
        axis.scatter(gdat.anglfact * lgal[indx], gdat.anglfact * bgal[indx], s=mrkrsize[indx], alpha=gdat.mrkralph, label=gdat.truelablhits, marker='*', linewidth=2, color='g')
        
        ## annotate
        # temp
        gdat.boolanot = False
        if gdat.boolanot:
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
    mrkrsize = retr_mrkrsize(gdat, gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[indxpoplplot][gdat.indxenerfluxdist, :]])
    lgal = gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[indxpoplplot]]
    bgal = gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[indxpoplplot]]
    axis.scatter(gdat.anglfact * lgal, gdat.anglfact * bgal, s=mrkrsize, alpha=gdat.mrkralph, label='Sample', marker='+', linewidth=2, color='b')


def retr_levi(listllik):
    
    minmlistllik = amin(listllik)
    levi = log(mean(1. / exp(listllik - minmlistllik))) + minmlistllik
    
    return levi


def retr_info(listllik, levi):
    
    info = mean(listllik) - levi

    return info


def retr_imag(gdat, axis, maps, thisindxener, thisindxevtt, logt=False, cmap='Reds', mean=False, satuuppr=None, satulowr=None, titl=''):

    if logt:
        maps = log10(abs(copy(maps)) + 1e-10)

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
    
    # saturate the map
    if satulowr != None:
        maps[where(maps < satulowr[thisindxener])] = satulowr[thisindxener]
    if satuuppr != None:
        maps[where(maps > satuuppr[thisindxener])] = satuuppr[thisindxener]
   
    # temp
    #imag = axis.imshow(maps)
    imag = axis.imshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='nearest')
    #imag = axis.matshow(maps, cmap=cmap, origin='lower', extent=gdat.exttrofi, interpolation='none')
    axis.set_title(titl)

    # make a color bar
    if thisindxevtt != None or thisindxener != None:
        cbar = plt.colorbar(imag, ax=axis, fraction=0.05)

    return axis, cbar


def retr_jcbn():
   
    fluxinit, lgalinit, bgalinit, sindinit, fluxauxi, radiauxi, anglauxi, sindauxi \
                                                            = sympy.symbols('fluxinit lgalinit bgalinit sindinit fluxauxi radiauxi anglauxi sindauxi')
    
    matr = sympy.Matrix([[     fluxauxi, 0, 0 , 0,                         fluxinit,                                    0,                                               0, 0], \
                         [            0, 1, 0 , 0, -radiauxi * sympy.cos(anglauxi) , (1 - fluxauxi) * sympy.cos(anglauxi), (fluxauxi - 1) * radiauxi * sympy.sin(anglauxi), 0], \
                         [            0, 0, 1 , 0, -radiauxi * sympy.cos(anglauxi) , (1 - fluxauxi) * sympy.cos(anglauxi), (1 - fluxauxi) * radiauxi * sympy.cos(anglauxi), 0], \
                         [            0, 0, 0 , 1,                                0,                                    0,                                               0, 0], \
                         [ 1 - fluxauxi, 0, 0 , 0,                        -fluxinit,                                    0,                                               0, 0], \
                         [            0, 1, 0 , 0, -radiauxi * sympy.cos(anglauxi) , -fluxauxi * sympy.cos(anglauxi), fluxauxi * radiauxi * sympy.sin(anglauxi), 0], \
                         [            0, 0, 1 , 0, -radiauxi * sympy.cos(anglauxi) , -fluxauxi * sympy.cos(anglauxi), -fluxauxi * radiauxi * sympy.cos(anglauxi), 0], \
                         [            0, 0, 0 , 0,                               0 ,                                    0,                                               0, 1]])

    jcbn = matr.det()
    print jcbn

    return jcbn


def corr_catl(gdat, thisindxpopl, modllgal, modlbgal, modlspec):

    indxtruepntsassc = gdatstrt()
    indxtruepntsassc.miss = []
    indxtruepntsassc.bias = [[] for i in gdat.indxener]
    indxtruepntsassc.hits = [[] for i in gdat.indxener]
    indxtruepntsassc.mult = []
        
    # get the flux limit that delineates the biased associations and hits 
    fluxbias = empty((2, gdat.numbener, gdat.numbflux + 1))
    for i in gdat.indxener:
        fluxbias[:, i, :] = retr_fluxbias(gdat, gdat.binsspec, i)

    indxmodlpnts = zeros_like(gdat.truelgal[thisindxpopl], dtype=int) - 1
    specassc = zeros((gdat.numbener, gdat.truenumbpnts[thisindxpopl]), dtype=float)
    numbassc = zeros_like(gdat.truelgal[thisindxpopl], dtype=int)
    distassc = zeros_like(gdat.truelgal[thisindxpopl]) + 3 * gdat.maxmgang
    dir1 = array([gdat.truelgal[thisindxpopl], gdat.truebgal[thisindxpopl]])

    for k in range(modllgal.size):
        dir2 = array([modllgal[k], modlbgal[k]])
        dist = angdist(dir1, dir2, lonlat=True)
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

    # check whether the flux of the associated model point source matches well with the flux of the deterministic point source
    for k in range(gdat.truenumbpnts[thisindxpopl]):
        if numbassc[k] == 0:
            indxtruepntsassc.miss.append(k)
        else:
            if numbassc[k] > 1:
                indxtruepntsassc.mult.append(k)
            for i in gdat.indxener:
                fluxbiasthis = interp(gdat.truespec[thisindxpopl][0, i, k], gdat.binsspec[i, :], fluxbias[0, i, :])
                boolbias = specassc[i, k] > fluxbiasthis or specassc[i, k] < gdat.truespec[thisindxpopl][0, i, k]**2 / fluxbiasthis 
                if boolbias:
                    indxtruepntsassc.bias[i].append(k)
                else:
                    indxtruepntsassc.hits[i].append(k)
   
    return indxmodlpnts, indxtruepntsassc


def retr_fluxbias(gdat, spec, indxenerthis):

    fluxbias = empty((2, spec.shape[1]))
    deno = log10(gdat.maxmspec[indxenerthis]) - log10(gdat.minmspec[indxenerthis])
    
    factlowr = 10.
    factuppr = 1.1
    
    offs = (log10(gdat.maxmspec[indxenerthis]) * log10(factlowr) - log10(gdat.minmspec[indxenerthis]) * log10(factuppr)) / deno
    
    slop = (log10(factuppr) + log10(gdat.maxmspec[indxenerthis]) - log10(factlowr) - log10(gdat.minmspec[indxenerthis])) / deno
    fluxbias[0, :] = 10**(slop * log10(spec[indxenerthis, :]) + offs)

    slop = (-log10(factuppr) + log10(gdat.maxmspec[indxenerthis]) + log10(factlowr) - log10(gdat.minmspec[indxenerthis])) / deno
    fluxbias[1, :] = 10**(slop * log10(spec[indxenerthis, :]) - offs)

    return fluxbias


