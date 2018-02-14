from __init__ import *
from util import *

print 'Bootstrapping PCAT runs over tiles...'

strgcnfg = sys.argv[1]
pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'

print 'Filter:'
print strgcnfg

listrtag = fnmatch.filter(os.listdir(pathdata), strgcnfg)

if '20180207_002838_pcat_chan_inpt_home7msc06000024_150000' in listrtag:
    print 'Removing 20180207_002838_pcat_chan_inpt_home7msc06000024_150000...'
    listrtag.remove('20180207_002838_pcat_chan_inpt_home7msc06000024_150000')

if len(listrtag) == 0:
    print 'Did not find any run tags.'
else:
    print 'Found the following run tags: '
    for rtag in listrtag:
        print rtag
    proc_finl(rtag=listrtag)


