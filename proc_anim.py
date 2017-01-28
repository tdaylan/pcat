import os
path = '/Users/tansu/Documents/work/git/local/docs/pres/20170128_dcol'
cntr = 0
for strg in os.listdir(path):
    if strg.startswith('thismodlcnts') and strg.endswith('pdf'):
        strgtemp = strg[:-4]
        cmnd = 'convert %s.pdf thismodlcnts%d.png' % (strgtemp, cntr)
        print cmnd
        os.system(cmnd)
        cntr += 1

