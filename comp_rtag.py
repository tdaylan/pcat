from __init__ import *

def comp(nameplot):
    
    if nameplot == 'DS_Store':
        return
    
    nameplot = nameplot[:-4]

    print 'comp() working on ' + nameplot + '...'
    
    cmnd = 'mkdir -p ' + pathimag + 'comprtag/' + nameplot + '/'
    print cmnd
    os.system(cmnd)
    
    print 'listline'
    print listline
    for line in listline:
        
        strgtemp = 'rtag'
        namedest = pathimag + 'comp' + strgtemp + '/' + nameplot + '/' + line + '.pdf'
        if not os.path.isfile(namedest):
            if boolfink:
                cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink2/www/tansu/link/pcat/imag/' + line + '/' + nameplot + '.pdf ' + namedest 
            else:
                cmnd = 'cp %s/' % pathimag + line + '/' + nameplot + '.pdf ' + namedest 
            
            print cmnd
            print
            os.system(cmnd)


print 'comp_rtag() initialized...'

if os.uname()[1] == 'fink1.rc.fas.harvard.edu' or os.uname()[1] == 'fink2.rc.fas.harvard.edu':
    boolfink = True
else:
    boolfink = False

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'

print 'Listing the available runs...'
if boolfink:
    cmnd = 'ssh tansu@fink2.rc.fas.harvard.edu "ls /n/fink1/tansu/data/pcat/imag/ | xargs -n 1 basename > /n/fink1/tansu/data/pcat/imag/listfink.txt"'
else:
    cmnd = 'ls %s | xargs -n 1 basename > %s/listlocl.txt' % (pathimag, pathimag)
print cmnd
os.system(cmnd)

if boolfink:
    print 'Copying the list to the local directory...'
    cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/imag/listfink.txt %s' % pathimag
    print cmnd
    os.system(cmnd)

if boolfink:
    pathlist = pathimag + 'listfink.txt'
else:
    pathlist = pathimag + 'listlocl.txt'

with open(pathlist) as thisfile:
    listline = thisfile.readlines()
    listline = [x.strip() for x in listline] 
strgsrch = '20180129_*_test_info_*_100000'
listline = fnmatch.filter(listline, strgsrch)

print 'listline'
for line in listline:
    print line

namerefr = '20180125_222826_pcat_chan_inpt_home7msc06000000_10000'
pathrefr = os.environ["PCAT_DATA_PATH"] + '/imag/' + namerefr + '/'

if boolfink:
    if not os.path.isdir(pathrefr):
        print 'Copying reference run to the local directory...'
        cmnd = 'scp -r tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/imag/%s/ %s' % (namerefr, pathimag)
        print cmnd
        os.system(cmnd)
        print

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



