from __init__ import *
from util import *

rtag = sys.argv[1]
strgplot = sys.argv[2]
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/%s/post/fram/' % rtag
listfram = fnmatch.filter(os.listdir(pathimag), '*%s*.pdf' % strgplot)
cmnd = 'convert '
numbfram = len(listfram)
indxfram = arange(numbfram)
indxframrand = choice(indxfram, size=numbfram, replace=False)
for fram in listpathfram[indxframrand]:
    cmnd += '%s ' % fram
cmnd += '%s.gif' % strgplot
#os.system(cmnd)
print cmnd

