import pcat.main, os, numpy as np, pyfits as pf

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
                   truelablback=['isottuto.fits'], \
                   )
                
pcat_tuto()

