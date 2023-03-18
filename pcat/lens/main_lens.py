

def plot_lens(gdat):
    
    if gmod.boolelemdeflsubh:
        xdat = gdat.binspara.angl[1:] * gdat.anglfact
        lablxdat = gdat.labltotlpara.gang
        
        listdeflscal = np.array([4e-2, 4e-2, 4e-2]) / gdat.anglfact
        listanglscal = np.array([0.05, 0.1, 0.05]) / gdat.anglfact
        listanglcutf = np.array([1.,    1.,  10.]) / gdat.anglfact
        listasym = [False, False, False]
        listydat = []
        for deflscal, anglscal, anglcutf, asym in zip(listdeflscal, listanglscal, listanglcutf, listasym):
            listydat.append(chalcedon.retr_deflcutf(gdat.binspara.angl[1:], deflscal, anglscal, anglcutf, asym=asym) * gdat.anglfact)
        
        for scalxdat in ['self', 'logt']:
            path = gdat.pathinitintr + 'deflcutf' + scalxdat + '.%s' % gdat.typefileplot
            tdpy.plot_gene(path, xdat, listydat, scalxdat=scalxdat, scalydat='logt', lablxdat=lablxdat, \
                                                             lablydat=r'$\alpha_n$ [$^{\prime\prime}$]', limtydat=[1e-3, 1.5e-2], limtxdat=[None, 2.])
       
        # pixel-convoltuion of the Sersic profile
        # temp -- y axis labels are wrong, should be per solid angle
        xdat = gdat.binspara.xpossers * gdat.anglfact
        for n in range(gdat.numbindxsers + 1):
            for k in range(gdat.numbhalfsers + 1):
                if k != 5:
                    continue
                path = gdat.pathinitintr + 'sersprofconv%04d%04d.%s' % (n, k, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, gdat.sersprof[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                #path = gdat.pathinitintr + 'sersprofcntr%04d%04d.%s' % (n, k, gdat.typefileplot)
                #tdpy.plot_gene(path, xdat, gdat.sersprofcntr[:, n, k], scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e6, 1e12])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.%s' % (n, k, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
                path = gdat.pathinitintr + 'sersprofdiff%04d%04d.%s' % (n, k, gdat.typefileplot)
                tdpy.plot_gene(path, xdat, abs(gdat.sersprof[:, n, k] - gdat.sersprofcntr[:, n, k]) / gdat.sersprofcntr[:, n, k], scalxdat='logt', \
                                                                     scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, limtydat=[1e-6, 1.])
       
        xdat = gdat.binspara.angl * gdat.anglfact
        listspec = np.array([1e-19, 1e-18, 1e-18, 1e-18]) / gdat.anglfact
        listsize = np.array([0.3, 1., 1., 1.]) / gdat.anglfact
        listindx = np.array([4., 2., 4., 10.])
        listydat = []
        listlabl = []
        for spec, size, indx in zip(listspec, listsize, listindx):
            listydat.append(spec * retr_sbrtsersnorm(gdat.binspara.angl, size, indxsers=indx))
            listlabl.append('$R_e = %.3g ^{\prime\prime}, n = %.2g$' % (size * gdat.anglfact, indx))
        path = gdat.pathinitintr + 'sersprof.%s' % gdat.typefileplot
        tdpy.plot_gene(path, xdat, listydat, scalxdat='logt', scalydat='logt', lablxdat=lablxdat, lablydat=gdat.lablfluxtotl, \
                                                                                                   listlegd=listlegd, listhlin=1e-7, limtydat=[1e-8, 1e0])
    
        #valulevl = np.linspace(7.5, 9., 5)
        valulevl = [7.0, 7.3, 7.7, 8., 8.6]
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        cont = axis.contour(gdat.binspara.redshost, gdat.binspara.redssour, minmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=20, fmt='%.3g')
        axis.set_xlabel(r'$z_{\rm{hst}}$')
        axis.set_ylabel(r'$z_{\rm{src}}$')
        axis.set_title(r'$M_{c,min}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsminm.%s' % gdat.typefileplot
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
        
        valulevl = np.linspace(9., 11., 20)
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        imag = axis.imshow(maxmmass, extent=[minmredshost, maxmredshost, minmredssour, maxmredssour], aspect='auto', vmin=9., vmax=11.)
        cont = axis.contour(gdat.binspara.redshost, gdat.binspara.redssour, maxmmass, 10, colors='g', levels=valulevl)
        axis.clabel(cont, inline=1, fontsize=15, fmt='%.3g')
        axis.set_xlabel('$z_{hst}$')
        axis.set_ylabel('$z_{src}$')
        axis.set_title(r'$M_{c,max}$ [$M_{\odot}$]')
        path = gdat.pathinitintr + 'massredsmaxm.%s' % gdat.typefileplot
        plt.colorbar(imag) 
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        axis.plot(gdat.meanpara.redshost, gdat.meanpara.adislens * gdat.sizepixl * 1e-3)
        axis.plot(gdat.meanpara.redshost, gdat.meanpara.adislens * 2. * gdat.maxmgangdata * 1e-3)
        axis.set_xlabel('$z_h$')
        axis.set_yscale('log')
        axis.set_ylabel(r'$\lambda$ [kpc]')
        path = gdat.pathinitintr + 'wlenreds.%s' % gdat.typefileplot
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots(figsize=(gdat.plotsize, gdat.plotsize))
        fracacutasca = np.logspace(-1., 2., 20)
        mcut = retr_mcutfrommscl(fracacutasca)
        axis.lognp.log(fracacutasca, mcut)
        axis.set_xlabel(r'$\tau_n$')
        axis.set_ylabel(r'$M_{c,n} / M_{0,n}$')
        axis.axhline(1., ls='--')
        path = gdat.pathinitintr + 'mcut.%s' % gdat.typefileplot
        plt.tight_layout()
        print('Writing to %s...' % path)
        figr.savefig(path)
        plt.close(figr)
       


