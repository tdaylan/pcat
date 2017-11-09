from __init__ import *

pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'

cmnd = 'rm -rf %s%s' % (pathdata, sys.argv[1])
print cmnd
os.system(cmnd)
cmnd = 'rm -rf %s%s' % (pathimag, sys.argv[1])
os.system(cmnd)
print cmnd
