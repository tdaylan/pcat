import os, sys
print 'sys.argv'
print sys.argv
strgtotl = ''
for strg in sys.argv[1:]:
    strgtotl += ' ' + strg
cmnd = 'python ' + strgtotl + ' > /dev/null 2>&1 &'

print cmnd
os.system(cmnd)
