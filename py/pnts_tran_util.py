
# coding: utf-8

# In[13]:

# math
from numpy import *
from scipy import interpolate

# astro & healpix
import healpy as hp
import pyfits as pf

# utilities
import glob
import os

# plotting
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='dark', context='poster', color_codes=True)

# tdpy
import tdpy_util


# In[18]:

def test_data():
    
    nside = 256
    npixl = nside**2 * 12
    evtt = [4,8,16,32]
    evtc = 512

    # energy axis for comparison with the older maps
    #ener, binsener, diffener, nener, nevtt = pnts_tran_util.get_axes(medi=True)
    #enerfull, binsenerfull, diffenerfull, nenerfull, nevtt = pnts_tran_util.get_axes(full=True)

    # new maps binned using pnts_tran energy axis
    expo = zeros((nener, npixl, nevtt))
    flux = zeros((nener, npixl, nevtt))    
    for m in range(nevtt):
        filelist = flist("medi", nside, evtc, evtt[m])
        for i in range(nener):
            expo[i,:,m] = pf.getdata(filelist[i], 1) # [cm^2*s] 
            flux[i,:,m] = pf.getdata(filelist[i], 0) / diffener[i] * 1e4 # [1/cm^2/s/sr/GeV]

            
    # P8R2_V6 maps in the old energy bins for comparison
    expofull = zeros((nenerfull, npixl, nevtt))
    fluxfull = zeros((nenerfull, npixl, nevtt))
    for m in range(nevtt):
        filelist = flist("full", nside, evtc, evtt[m])
        for i in range(nenerfull):
            expofull[i,:,m] = pf.getdata(filelist[i], 1) # [cm^2*s] 
            fluxfull[i,:,m] = pf.getdata(filelist[i], 0) / diffenerfull[i] * 1e4 # [1/cm^2/s/sr/GeV]
    expofulltotl = sum(expofull, axis=2)
    fluxfullmean = sum(fluxfull * expofull, axis=2) / expofulltotl
    
    
    # P7REP_V15 maps
    evttrep7 = ['p17_ultraclean', 'p17_ultraclean_Q2', 'p17_ultraclean_Q1']
    nevttrep7 = len(evttrep7)
    exporep7 = zeros((nenerfull, npixl, nevttrep7))
    fluxrep7 = zeros((nenerfull, npixl, nevttrep7))
    for k in range(nevttrep7):
        for i in range(nenerfull):
            filelist = flist_legacy(evttrep7[k], nopsc=True)
            exporep7[i,:,k] = pf.getdata(filelist[i], 1) # [cm^2*s] 
            fluxrep7[i,:,k] = pf.getdata(filelist[i], 0) / diffenerfull[i] * 1e4 # [1/cm^2/s/sr/GeV]
    

    lonra = [-20., 20.]
    latra = [-10., 10.]
    extent = [lonra[1], lonra[0], latra[0], latra[1]]
    fig, axgrd = plt.subplots(8, 3, figsize=(20, 30))
    titles = ['P8R2_V6', 'P7REP_V15', 'Residual']
    fig.suptitle(r'Flux [10$^4$ $\times$ GeV/cm$^2$/s/sr], %.3g GeV - %.3g GeV' % (binsenerfull[0], binsenerfull[1]), fontsize=25)
    for r, axrow in enumerate(axgrd):
        for n, ax in enumerate(axrow):
            if r < 4:
                if n == 0:
                    title = 'P8R2_V6 PSF%d' % r
                    data = fluxfull[0,:,r]
                if n == 1:
                    title = 'P7REP_V15'
                    data = fluxrep7[0,:,0]
                if n == 2:
                    title = 'P8R2_V6 PSF%d - P7REP_V15' % r
                    data = fluxfull[0,:,r] - fluxrep7[0,:,0]
                    
            if r == 4:
                if n == 0:
                    title = 'P8R2_V6 AVG'
                    data = fluxfullmean[0,:]
                if n == 1:
                    title = 'P7REP_V15'
                    data = fluxrep7[0,:,0]
                if n == 2:
                    title = 'P8R2_V6 AVG - P7REP_V15'
                    data = fluxfullmean[0,:] - fluxrep7[0,:,0]
                    
            if r == 5:
                if n == 0:
                    title = 'P8R2_V6 AVG'
                    data = hp.smoothing(fluxfullmean[0,:], pi / 90., verbose=False)
                if n == 1:
                    title = 'P7REP_V15'
                    data = hp.smoothing(fluxrep7[0,:,0], pi / 90., verbose=False)
                if n == 2:
                    title = 'P8R2_V6 AVG - P7REP_V15'
                    data = hp.smoothing(fluxfullmean[0,:] - fluxrep7[0,:,0], pi / 90., verbose=False)
                
                
            if r == 6:
                if n == 0:
                    title = r'P8R2_V6 Q$^\prime$1'
                    data = fluxfull[0,:,3]
                if n == 1:
                    title = 'P7REP_V15 Q1'
                    data = fluxrep7[0,:,2]
                if n == 2:
                    title = r'P8R2_V6 Q$^\prime$1 - P7REP_V15 Q1'
                    data = fluxfull[0,:,3] - fluxrep7[0,:,2]
                    
            if r == 7:
                if n == 0:
                    title = r'P8R2_V6 Q$^\prime$2'
                    data = (fluxfull[0,:,3] * expofull[0,:,3] + fluxfull[0,:,2] * expofull[0,:,2]) / (expofull[0,:,3] + expofull[0,:,2])
                if n == 1:
                    title = 'P7REP_V15 Q2'
                    data = fluxrep7[0,:,1]
                if n == 2:
                    title = r'P8R2_V6 Q$^\prime$2 - P7REP_V15 Q2'
                    data = (fluxfull[0,:,3] * expofull[0,:,3] + fluxfull[0,:,2] * expofull[0,:,2]) / (expofull[0,:,3] + expofull[0,:,2]) - fluxrep7[0,:,1]
                    
            ax.set_title(title)
            im = tdpy_util.get_cart(data, latra=latra, lonra=lonra)
            im = ax.imshow(im, cmap='Reds', extent=extent)
            plt.colorbar(im, ax=ax, fraction=0.02)
    fig.tight_layout()
    fig.subplots_adjust(top=0.95, hspace=0.2)
    plt.show()   
    
    expo, exprflux = get_expoflux(256)
        
    lonra = [-180, 180]
    latra = [-90, 90]
    extent = lonra[1], lonra[0], latra[0], latra[1]
    
    fig, axgrd = plt.subplots(nevtt, nener, figsize=(20, 20))
    fig.suptitle('Exposure [cm$^2$ s]', fontsize=25)
    for m, axrow in enumerate(axgrd):
        for i, ax in enumerate(axrow):
            testcart = tdpy_util.get_cart(expo[i,:,m], latra=latra, lonra=lonra)
            imag = ax.imshow(testcart, origin='lower', cmap='Reds', extent=extent)
            plt.colorbar(imag, ax=ax, fraction=0.04)
            if m == 0:
                ax.set_title('%.3g GeV - %.3g GeV' % (binsener[i], binsener[i]))
            if i == 0:
                ax.set_ylabel('PSF %d' % m)
    fig.subplots_adjust(top=0.95, hspace=0.25)
    plt.show()


# In[20]:

def flist(version, nside, evclass, evtype):

    tag = os.environ["PNTS_TRAN_DATA_PATH"] + "/" + version + "/flux_*_%4.4d_evc%3.3d_evt%2.2d.fits" % (nside, evclass, evtype)
    flist = sorted(glob.glob(tag))

    return flist

