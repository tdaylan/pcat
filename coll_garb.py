from __init__ import *

pathchec = os.environ["PCAT_DATA_PATH"] + '/data/outp/comp.txt'

if len(sys.argv) == 2 and sys.argv[1] == 'data':
    liststrgextn = ['/data/outp/']
else:
    liststrgextn = ['/imag/', '/data/outp/']

for strgextn in liststrgextn:

    path = os.environ["PCAT_DATA_PATH"] + strgextn

    for strgfile in os.listdir(path):
        pathfile = path + strgfile
        if os.path.isdir(pathfile) and strgfile[:8].isdigit():

            pathchec = pathfile.replace('imag', 'data/outp') + '/comp.txt'

            try:
                numbswep = int(pathfile[pathfile.rfind('_')+1:])
                if not os.path.exists(pathchec) or numbswep < 1000000:
                    os.system('rm -rf ' + pathfile)
                    print 'Deleting %s' % pathchec
                else:
                    print 'Saving %s' % pathchec
            except:
                pass

