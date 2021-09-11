import os
import sys

import numpy as np

import pcat

def pcat_tuto():
    # you should set $PCAT_DATA_PATH to a folder that will contain input, indermediate and output files of PCAT
    #os.environ['PCAT_DATA_PATH'] = '/path/to/your/pcat/data/folder'
    
    numbener = 5
    numbside = 256
    numbpixl = 12 * numbside**2
    numbevtt = 4
    fluxisot = 1e-6 * np.ones((numbener, numbpixl, numbevtt))
    
    # write the isotropic template to the PCAT input folder
    # Note that the default expected name is fermisotflux.fits
    # but could be different if strgback argument was provided 
    path = os.environ['PCAT_DATA_PATH'] + '/data/inpt/isottuto.fits'
    pf.writeto(path, fluxisot, clobber=True)
    
    # For this run, we will assume a Gaussian prior on the PSF parameters, centered on the best-fit Fermi-LAT PSF 
    # Therefore, we need to download the Fermi-LAT IRF files.
    cmnd = 'wget https://faun.rc.fas.harvard.edu/tansu/pcat/tuto/psf_P7REP_SOURCE_V15_back.fits $PCAT_DATA_PATH/data/inpt/psf_P7REP_SOURCE_V15_back.fits'
    #os.system(cmnd)
    
    # Now we can run PCAT
    pcat.init( \
                   forccart=True, \
                   pixltype='cart', \
                   diagmode=False, \
                   backtype=[1.], \
                   numbswep=2000000, \
                   strgexpo=1e11, \
                   probbrde=0.5, \
                   #lablback=['Isotropic'], \
                   )
                

def retr_modl(para, grid):
    
    parafixp, paraelem = para
    
    modl = np.sum(paraelem['ampl'][None, :] * np.exp(-0.5 * (paraelem['posi'][None, :] - grid[:, None])**2), 1)
    
    return modl


def retr_llik(para, *args):
    
    modl = retr_modl(para, grid)

    llik = np.sum(-0.5 * (data - modl) / vari)
    
    return llik


def cnfg_GMM():

    # generate grid
    numbgrid = 100
    grid = np.linspace(-1, 1., numbgrid)
    
    # make mock data
    numbelem = 10
    parafixp = np.array([1.])
    paraelem = {}
    paraelem['ampl'] = np.random.rand(numbelem)
    paraelem['posi'] = np.random.rand(numbelem) * 2. - 1.
    para = parafixp, paraelem
    obsd = retr_modl(para, grid) + 1e-2 * np.random.randn(numbgrid)
    
    pcat.init( \
              typeexpr='gene', \
              typeverb=0, \
             )


def cnfg_external():

    # generate synethetic data
    ## number of samples
    numbsamp = 100
    ## horizontal (x-axis) position of the synethetic samples
    xpossamp = pcat.icdf_gaus(np.random.rand(numbsamp), 0.1, 0.5)
    ## vertical (y-axis) position of the synethetic samples
    ypossamp = pcat.icdf_gaus(np.random.rand(numbsamp), 0.2, 0.1)

    def retr_llik(para, paraelem, *args):
        
        xpossamp, ypossamp = args

        llik = np.exp(-(xpos - xpossamp)**2) + np.exp(-(yposgrid - ypossamp)**2)
        
        return llik


    def retr_lpri(para, paraelem, *args):
        
        lpri = np.exp(-xpos**2) + np.exp(-yposgrid**2)
        
        return lpri

    pcat.init(retr_llik=retr_llik, retr_lpri=retr_lpri)


globals().get(sys.argv[1])()

