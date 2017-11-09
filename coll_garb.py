from __init__ import *

pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/'
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'

listrtagdata = fnmatch.filter(os.listdir(pathdata), '2*')
listrtagimag = fnmatch.filter(os.listdir(pathimag), '2*')

for rtag in listrtagdata:
    if not rtag in listrtagimag:
        raise Exception('Data and image folders are not synched!')
for rtag in listrtagimag:
    if not rtag in listrtagdata:
        raise Exception('Data and image folders are not synched!')

liststrgextn = ['/imag/', '/data/outp/']

for strgextn in liststrgextn:

    path = os.environ["PCAT_DATA_PATH"] + strgextn
    
    for rtag in os.listdir(path):
        

        pathfile = path + rtag
        if os.path.isdir(pathfile) and rtag[:8].isdigit():
            print 'Processing %s...' % rtag

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
                if (not os.path.isfile(pathchec) or not booltemp or numbswep < 1000) and not 'mockonly' in rtag:
                    if not rtag == '20171025_115229_pcat_lens_mock_syst_lowrtrue_2000000':
                        cmnd = 'rm -rf ' + pathfile
                        os.system(cmnd)
                    else:
                        print 'Saving %s' % rtag
                    print 'Deleting %s' % pathchec
                else:
                    pass
                    print 'Saving %s' % pathchec
            except:
                print 'Skipping...'
                pass
            print

