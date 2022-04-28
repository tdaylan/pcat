import os
import sys

import numpy as np

import tdpy
import pcat

def cnfg_fermi():
    '''
    The output folder will contain initial, indermediate, and output image and data files. To define this folder, tou should either
        1) set $PCAT_DATA_PATH environment variable, as in
            os.environ['PCAT_DATA_PATH'] = '/path/to/your/pcat/data/folder'
            inside your .bashrc or .zshrc, or
        2) input pathpcat as an argument.
    '''
    
    
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
                

def cnfg_GaussianMix():
    '''
    Gaussian Mixture
    '''
    
    pcat.init( \
             )


def cnfg_GaussianMix_unbinned():
    '''
    Unbinned Gaussian Mixture
    '''
    
    pcat.init( \
              boolbins=False, \
             )


globals().get(sys.argv[1])()

