from __init__ import *
from util import *

rtagroot = sys.argv[1]
pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'

listrtagdata = fnmatch.filter(os.listdir(pathdata), rtagroot)
numbcnfg = len(listrtagdata)
for k, rtag in enumerate(listrtagdata):
    print 'Processing %s...' % rtag
    pathoutprtag = retr_pathoutprtag(rtag)
    path = pathoutprtag + 'gdatinit'
    gdat = readfile(path) 
    if k == 0:
        cntpdataarry = empty([numbcnfg] + list(gdat.cntpdatareg0.shape))
    cntpdataarry[k, ...] = gdat.cntpdatareg0

path = pathdata + 'truecntpdata_%s.h5' % listrtagdata[0]
print 'Writing to %s...' % path
filearry = h5py.File(path, 'w')
filearry.create_dataset('cntpdataarry', data=cntpdataarry)
filearry.close()


