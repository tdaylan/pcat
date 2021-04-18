from __init__ import *
from util import *

print 'PCAT plot merging routine'

if len(sys.argv) > 2:
    listrtag = sys.argv[1:]
else:
    strgsrch = sys.argv[1]
    pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
    listrtag = os.listdir(pathimag)
    listrtag = fnmatch.filter(listrtag, strgsrch)

listgdat = retr_listgdat(listrtag, typegdat='init')

for gdat in listgdat:
    merg_plot(gdat)
    
