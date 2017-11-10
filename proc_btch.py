from __init__ import *
from util import *

print 'PCAT post process routine started.'

boolfinl = False
boolanim = False
if len(sys.argv) == 1:
    strgsrch = '20*'
    boolfinl = True
    boolanim = True
elif len(sys.argv) == 2:
    strgsrch = sys.argv[1]
    boolfinl = True
    boolanim = True
elif len(sys.argv) == 3:
    strgsrch = sys.argv[1]
    boolfinl = sys.argv[2] == 'finl'
    boolanim = sys.argv[2] == 'anim'
    
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listrtag = fnmatch.filter(os.listdir(pathimag), strgsrch)

if boolfinl:
    print 'Post-processing...'
    for rtag in listrtag:
        try:
            print 'Working on %s...' % rtag
            proc_finl(rtag=rtag)
        except:
            pass
    
if boolanim:
    print 'Making animations...'
    for rtag in listrtag:
        print 'Working on %s...' % rtag
        proc_anim(rtag=rtag)
    
