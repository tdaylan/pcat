import os, sys
from __init__ import * 

strgtotl = ''
for strg in sys.argv[1:]:
    strgtotl += ' ' + strg

strgtimestmp = tdpy.util.retr_strgtimestmp()

cmnd = 'python' + strgtotl + ' > ' + os.environ["PCAT_DATA_PATH"] + '/data/' + strgtimestmp + '_' + sys.argv[2] + '_' + sys.argv[3] + '.out 2>&1 &'

print 'cmnd'
print cmnd

os.system(cmnd)
