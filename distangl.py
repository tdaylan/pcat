from numpy import *
from numpy.random import *
import time
import healpy as hp

def retr_unit(lgal, bgal):

    lgaltemp = deg2rad(lgal)
    bgaltemp = deg2rad(bgal)

    xaxi = cos(bgaltemp) * cos(lgaltemp)
    yaxi = cos(bgaltemp) * sin(lgaltemp)
    zaxi = sin(bgaltemp)

    return xaxi, yaxi, zaxi


def retr_angldistunit(lgal, bgal):
    
    xaxi, yaxi, zaxi = retr_unit(lgal, bgal)
    acostemp = xaxigrid[None, :] * xaxi[:, None] + yaxigrid[None, :] * yaxi[:, None] + zaxigrid[None, :] * zaxi[:, None]
    
    angl = arccos(acostemp)

    return angl


def retr_angldistheal(lgal0, bgal0, lgal1, bgal1, aprx=False):
    
    if not aprx:
        dir1 = array([lgal0, bgal0])
        dir2 = array([lgal1, bgal1])
        dist = hp.rotator.angdist(dir1, dir2, lonlat=True)
    else:
        dist = deg2rad(sqrt((lgal0 - lgal1)**2 + (bgal0 - bgal1)**2))

    return dist

numbpnts = 2
numbpixl = 4
lgalgrid = linspace(-20., 20, numbpixl)
bgalgrid = linspace(-20., 20, numbpixl)

xaxigrid, yaxigrid, zaxigrid = retr_unit(lgalgrid, bgalgrid)

lgal = rand(numbpnts) * 40. - 20.
bgal = rand(numbpnts) * 40. - 20.

angl = empty(numbpnts)

numbtime = 10
listtime = empty((numbtime, 3))

for t in range(numbtime):
    
    timeinit = time.time()
    angl = retr_angldistunit(lgal, bgal)
    listtime[t, 0] = time.time() - timeinit
    print angl
    
    timeinit = time.time()
    for n in range(numbpnts):
        angl[n] = retr_angldistheal(lgal[n], bgal[n], lgalgrid, bgalgrid)
    listtime[t, 1] = time.time() - timeinit

    timeinit = time.time()
    for n in range(numbpnts):
        angl[n] = retr_angldistheal(lgal[n], bgal[n], lgalgrid, bgalgrid, aprx=True)
    listtime[t, 2] = time.time() - timeinit

print 'time'
print mean(listtime, 0) * 1e3




    

