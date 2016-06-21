import pyfits as pf, os
from numpy import *

pathsporcatl = os.environ["PCAT_DATA_PATH"] + '/sporcatl.fits'
listlgal = pf.getdata(pathsporcatl, 1).T
listbgal = pf.getdata(pathsporcatl, 2).T
listspec = pf.getdata(pathsporcatl, 3).T / array([0.7, 2., 7.])[None, :, None]

print listlgal[:, 0:4]
print listbgal[:, 0:4]
print listspec[:, 0, 0:4]
print listspec[:, 0, 0:4]
print listspec[:, 0, 0:4]

