from __init__ import *

def comp(nameplot):
    
    if nameplot == 'DS_Store':
        return
    
    nameplot = nameplot[:-4]

    print
    print 'comp() working on ' + nameplot + '...'
    
    cmnd = 'mkdir -p ' + pathimag + 'comprtag/' + nameplot + '/'
    print cmnd
    os.system(cmnd)
    
    for line in listline:
        
        strgtemp = 'rtag'
        namedest = pathimag + 'comp' + strgtemp + '/' + nameplot + '/' + line + '.pdf'
        if not os.path.isfile(namedest):
            cmnd = 'cp %s/' % pathimag + line + '/' + nameplot + '.pdf ' + namedest 
            
            print cmnd
            print
            os.system(cmnd)


print 'comp_rtag() initialized...'
print 'listline'
print listline

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'

print 'Listing the available runs...'
cmnd = 'ls %s | xargs -n 1 basename > %s/listrtag.txt' % (pathimag, pathimag)
print cmnd
os.system(cmnd)

pathlist = pathimag + 'listrtag.txt'

with open(pathlist) as thisfile:
    listline = thisfile.readlines()
    listline = [x.strip() for x in listline] 
strgsrch = sys.argv[1]
listline = fnmatch.filter(listline, strgsrch)

print 'listline'
for line in listline:
    print line

namerefr = listline[0]
pathrefr = os.environ["PCAT_DATA_PATH"] + '/imag/' + namerefr + '/'

print 'Iterating over folders of the reference run...'
for namefrst in os.listdir(pathrefr):
    namefrstprim = os.path.join(pathrefr, namefrst)
    if os.path.isdir(namefrstprim):
        for nameseco in os.listdir(namefrstprim):
            if nameseco == 'fram':
                continue
            namesecoprim = os.path.join(namefrstprim, nameseco)
            if os.path.isdir(namesecoprim):
                for namethrd in os.listdir(namesecoprim):
                    namethrdprim = os.path.join(namesecoprim, namethrd)
                    if os.path.isdir(namethrdprim):
                        for namefrth in os.listdir(namethrdprim):
                            namefrthprim = os.path.join(namethrdprim, namefrth)
                            comp(namefrthprim.split(pathrefr)[-1])
                    else:
                        comp(namethrdprim.split(pathrefr)[-1])
            elif nameseco.endswith('pdf') or nameseco.endswith('gif'):
                comp(namesecoprim.split(pathrefr)[-1])
    elif namefrst.endswith('pdf') or namefrst.endswith('gif'):
        comp(namefrstprim.split(pathrefr)[-1])



