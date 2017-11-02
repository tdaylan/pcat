from __init__ import *

if len(sys.argv) == 2 and sys.argv[1] == 'data':
    liststrgextn = ['/data/outp/']
else:
    liststrgextn = ['/imag/', '/data/outp/']

for strgextn in liststrgextn:

    path = os.environ["PCAT_DATA_PATH"] + strgextn
    
    for strgfile in os.listdir(path):
        pathfile = path + strgfile
        if os.path.isdir(pathfile) and strgfile[:8].isdigit():
            print 'Processing %s...' % strgfile

            # check the chain status
            pathchec = pathfile.replace('imag', 'data/outp') + '/stat.txt'
            if os.path.isfile(pathchec):
                filestat = open(pathchec, 'r')
                booltemp = False
                for line in filestat:
                    if line == 'gdatmodi written.\n':
                        booltemp = True
        
            try:
                numbswep = int(pathfile[pathfile.rfind('_')+1:])
                print 'numbswep'
                print numbswep
                print 'booltemp'
                print booltemp
                if not os.path.isfile(pathchec) or not booltemp or numbswep < 100000:
                    if not strgfile == '20171025_115229_pcat_lens_mock_syst_lowrtrue_2000000':
                        os.system('rm -rf ' + pathfile)
                    else:
                        print 'Saving %s' % strgfile
                    print 'Deleting %s' % pathchec
                else:
                    pass
                    print 'Saving %s' % pathchec
            except:
                print 'Skipping...'
                pass
            print

