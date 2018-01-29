from __init__ import *

def comp(nameplot):
    
    if nameplot == 'DS_Store':
        return
    
    nameplot = nameplot[:-4]

    print 'comp() working on ' + nameplot + '...'
    
    cmnd = 'mkdir -p ' + pathimag + 'comprtag/' + nameplot + '/'
    print cmnd
    #os.system(cmnd)
    
    cmnd = 'mkdir -p ' + pathimag + 'compgold/' + nameplot + '/'
    print cmnd
    #os.system(cmnd)
    
    for line in listline:
        
        if line in listlinegold:
            strgtemp = 'gold'
        else:
            strgtemp = 'rtag'
        namedest = pathimag + 'comp' + strgtemp + '/' + nameplot + '/' + line + '.pdf'
        if not os.path.isfile(namedest):
            if boolfink:
                cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink2/www/tansu/link/pcat/imag/' + line + '/' + nameplot + '.pdf ' + namedest 
            else:
                cmnd = 'cp %s/' % pathimag + line + '/' + nameplot + '.pdf ' + namedest 
            
            print cmnd
            #os.system(cmnd)


print 'comp_rtag() initialized...'

# list of golden runs
listlinegold = [ \
                '20180125_222826_pcat_chan_inpt_home7msc06000000_10000', \
                '20180125_223105_pcat_chan_inpt_home7msc06000001_10000', \
               ]

if os.uname()[1] == 'fink1.rc.fas.harvard.edu' or os.uname()[1] == 'fink2.rc.fas.harvard.edu':
    boolfink = True
else:
    boolfink = False

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
pathgold = pathimag + 'compgold/'

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
    #os.system(cmnd)

print 'Erasing the gold compilation contents...'
cmnd = 'rm -rf ' + pathgold
print cmnd
os.system(cmnd)

if boolfink:
    pathlist = pathimag + 'listfink.txt'
else:
    pathlist = pathimag + 'listlocl.txt'

with open(pathlist) as thisfile:
    listline = thisfile.readlines()
    listline = [x.strip() for x in listline] 
strgsrch = '20*chan_inpt*'
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
        #os.system(cmnd)

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



