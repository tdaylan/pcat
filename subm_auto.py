from __init__ import *

path = os.environ["TDGU_PATH"] + '/'
fileoutp = open(path + 'subm_auto.log', 'w')
cntr = 0
for name in os.listdir(path):
    if name.endswith(".py"):
        print name
        fileobjt = open(path + name, 'r')
        for line in fileobjt:
            if line.startswith('def pcat_'):
                
                #if cntr == 5:
                #    break
                
                cntr += 1
                namefunc = line[4:-1].split('(')[0]
                cmnd = 'python $TDGU_PATH/%s %s' % (name, namefunc)
                print cmnd
                try:
                    os.system(cmnd)
                    fileoutp.write('%s successfull.' % namefunc)
                except Exception as excp:
                    strg = str(excp)
                    fileoutp.write('%s failed.' % namefunc)
                    fileoutp.write(strg)
                print
                print
                print
                print
                print
                print
                print
                print

fileoutp.close()

