from __init__ import *

def comp(nameplot):
    
    if nameplot == 'DS_Store':
        return
    
    nameplot = nameplot[:-4]

    print 'comp() working on ' + nameplot + '...'
    
    cmnd = 'mkdir -p ' + pathimag + 'compfink/' + nameplot + '/'
    os.system(cmnd)
    
    cmnd = 'mkdir -p ' + pathimag + 'compgold/' + nameplot + '/'
    os.system(cmnd)
    
    for line in listline:
        
        if line in listlinegold:
            strgtemp = 'gold'
        else:
            strgtemp = 'fink'
        namedest = pathimag + 'comp' + strgtemp + '/' + nameplot + '/' + line + '.pdf'
        if not os.path.isfile(namedest):
            cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink2/www/tansu/link/pcat/imag/' + line + '/' + nameplot + '.pdf ' + namedest 
            os.system(cmnd)

print 'compfink initialized...'

# list of golden runs
listlinegold = [ \
                '20170522_020917_pcat_lens_mock_syst_0001_5000000', \
                '20170521_232045_pcat_lens_mock_syst_0000_5000000', \
               ]

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
pathgold = pathimag + 'compgold/'

print 'Listing the available runs...'
cmnd = 'ssh tansu@fink2.rc.fas.harvard.edu "ls /n/fink1/tansu/data/pcat/imag/ | xargs -n 1 basename > /n/fink1/tansu/data/pcat/imag/listfink.txt"'
os.system(cmnd)

print 'Copying the list to the local directory...'
cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/imag/listfink.txt %s' % pathimag
os.system(cmnd)

print 'Erasing the gold compilation contents...'
os.system('rm -rf ' + pathgold)

pathlist = pathimag + 'listfink.txt'
with open(pathlist) as thisfile:
    listline = thisfile.readlines()
    listline = [x.strip() for x in listline] 
strgsrch = '20170609*'
listline = fnmatch.filter(listline, strgsrch)

namerefr = '20170609_050550_pcat_lens_mock_syst_0000_2000000'
pathrefr = os.environ["PCAT_DATA_PATH"] + '/imag/' + namerefr + '/'
if not os.path.isdir(pathrefr):
    print 'Copying reference run to the local directory...'
    cmnd = 'scp -r tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/imag/%s/ %s' % (namerefr, pathimag)
    os.system(cmnd)

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



