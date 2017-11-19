from __init__ import *
from util import *

rtag = sys.argv[1]
strgplot = sys.argv[2]
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/%s/post/fram/' % rtag
listpathfram = array(fnmatch.filter(os.listdir(pathimag), '*%s*.pdf' % strgplot))
cmnd = 'convert '
numbfram = len(listpathfram)
indxfram = arange(numbfram)
indxframrand = choice(indxfram, size=numbfram, replace=False).astype(int)
for pathfram in listpathfram[indxframrand]:
    cmnd += '%s/%s ' % (pathimag, pathfram)
cmnd += '%s.gif' % strgplot
print cmnd
os.system(cmnd)

