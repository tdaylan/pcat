from numpy import *
import healpy as hp

fgl3bgal = deg2rad(5.)
fgl3lgal = deg2rad(5.)

rttr = hp.rotator.Rotator(rot=[fgl3bgal, fgl3lgal, 0.], deg=True)
fgl3bgal, fgl3lgal = rttr(pi / 2. - fgl3bgal, fgl3lgal)
fgl3bgal = pi / 2. - fgl3bgal

fgl3bgal *= 180. / pi
fgl3lgal *= 180. / pi


print 'fgl3lgal'
print fgl3lgal
print 'fgl3bgal'
print fgl3bgal
