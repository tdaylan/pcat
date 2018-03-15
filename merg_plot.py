from __init__ import *
from util import *

print 'PCAT plot merging routine'

if len(sys.argv) > 2:
    listrtag = sys.argv[1:]
else:
    strgsrch = sys.argv[1]
    listrtag = fnmatch.filter(listline, strgsrch)

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listrtag = fnmatch.filter(os.listdir(pathimag), '20*_pcat_*')
listgdat = retr_listgdat(listrtag, typegdat='init')

for gdat in listgdat:
    merg_plot(gdat)
    
