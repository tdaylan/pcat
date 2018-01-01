from __init__ import *
from util import *

print 'Bootstrapping PCAT runs over tiles...'

strgcnfg = sys.argv[1]
pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'

print 'Filter:'
print strgcnfg

listrtag = fnmatch.filter(os.listdir(pathdata), strgcnfg)

print 'Found the following run tags: '
for rtag in listrtag:
    print rtag

proc_finl(rtag=listrtag)


