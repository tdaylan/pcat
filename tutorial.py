import os

import numpy as np
import pyfits as pf

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
    pcat.main.init( \
                   forccart=True, \
                   pixltype='cart', \
                   diagmode=False, \
                   backtype=[1.], \
                   numbswep=2000000, \
                   strgexpo=1e11, \
                   probbrde=0.5, \
                   #lablback=['Isotropic'], \
                   )
                

def retr_modl(para):
    
    parafixp, paraelem = para
    
    modl = np.sum(paraelem['ampl'][None, :] * np.exp(-0.5 * (paraelem['posi'][None, :] - grid[:, None])**2), 1)
    
    return modl


def retr_llik(para, *args):
    
    modl = retr_modl(para)

    llik = np.sum(-0.5 * (data - modl) / vari)
    
    return llik


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
obsd = retr_modl(para) + 1e-2 * np.random.randn(numbgrid)

pcat.main.init(retr_llik=retr_llik)



