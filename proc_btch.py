from __init__ import *
from util import *

if sys.argv[1]:
    strgsrch = sys.argv[1]
else:
    strgsrch = '20*'
    
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listrtag = fnmatch.filter(os.listdir(pathimag), strgsrch)

print 'Batch production of animations...'

for rtag in listrtag:
    
    if sys.argv[2] == 'finl':
        proc_finl(rtag=rtag)
    
    if sys.argv[2] == 'anim':
        proc_anim(rtag=rtag)
    
