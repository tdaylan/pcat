from __init__ import *

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listpathruns = fnmatch.filter(os.listdir(pathimag), '20*')

print 'Making animations of frame plots...'

for pathruns in listpathruns:
    for namesampdist in ['prio', 'post']:
        for nameextn in ['', 'assc/', 'histodim/', 'histtdim/', 'scattdim/']:
            
            pathframextn = pathimag + pathruns + '/' + namesampdist + '/fram/' + nameextn
            pathanimextn = pathimag + pathruns + '/' + namesampdist + '/anim/' + nameextn
        
            try:
                listfile = fnmatch.filter(os.listdir(pathframextn), '*_swep*.pdf')
            except:
                print '%s failed.' % pathframextn
                continue

            listfiletemp = []
            for thisfile in listfile:
                listfiletemp.extend((thisfile.split('_')[0]).rsplit('/', 1))
            
            listname = list(set(listfiletemp))
            if len(listname) == 0:
                continue
            
            shuffle(listname)
            

            for name in listname:
                
                if not (name.startswith('thisdatacnts') or name.startswith('thishistdefspop0')):
                    continue

                strgtemp = '%s*_swep*.pdf' % name
                listfile = fnmatch.filter(os.listdir(pathframextn), strgtemp)
                numbfile = len(listfile)
                liststrgextn = []
                for k in range(numbfile):
                    liststrgextn.append((listfile[k].split(name)[1]).split('_')[0])
                
                liststrgextn = list(set(liststrgextn))
                
                for k in range(len(liststrgextn)):
            
                    listfile = fnmatch.filter(os.listdir(pathframextn), name + liststrgextn[k] + '_swep*.pdf')
                    numbfile = len(listfile)
                    
                    indxfilelowr = 0
                    
                    if indxfilelowr < numbfile:
                        indxfileanim = arange(indxfilelowr, numbfile)
                    else:
                        continue
                        
                    indxfileanim = choice(indxfileanim, replace=False, size=indxfileanim.size)
                    
                    cmnd = 'convert -delay 20 -density 300 -quality 100 '
                    for n in range(indxfileanim.size):
                        cmnd += '%s%s ' % (pathframextn, listfile[indxfileanim[n]])

                    namegiff = '%s%s.gif' % (pathanimextn, name + liststrgextn[k])
                    cmnd += ' ' + namegiff
                    if not os.path.exists(namegiff):
                        print 'Run: %s, pdf: %s' % (pathruns, namesampdist)
                        print 'Making %s animation...' % name
                        print
                        os.system(cmnd)
                    else:
                        pass

