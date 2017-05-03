from __init__ import *

path = '/Users/tansu/Downloads/pcatlite.fits'
tdpy.util.read_fits(path)

thisfile = pf.open(path)
lgal = thisfile[1].data
bgal = thisfile[2].data
flux = thisfile[3].data
sind = thisfile[4].data
meanpnts = thisfile[7].data
fluxdistslop = thisfile[8].data
psfp = thisfile[9].data
bacp = thisfile[10].data

set_printoptions(precision=2)
print lgal.shape
print lgal[:200, :20]
for k in range(lgal.shape[0]):
    print where(lgal[k, :] != 0.)[0].size
