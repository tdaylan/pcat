from __init__ import *
from util import *

strgcnfg = sys.argv[1]
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
print 'strgcnfg'
print strgcnfg
listrtag = fnmatch.filter(os.listdir(pathimag), strgcnfg)
print 'listrtag'
print listrtag
proc_finl(rtag=listrtag)


