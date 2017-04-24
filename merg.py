from __init__ import *

nameplot = 'gmrbmaps_00'

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listpathruns = fnmatch.filter(os.listdir(pathimag), '20*')
numbruns = len(listpathruns)
listnamesampdist = []
listpath = []
listrtag = []
for k in range(numbruns):
    for namesampdist in ['prio', 'post']:
        pathgmrb = pathimag + listpathruns[k] + '/' + namesampdist + '/finl/diag/' + nameplot + '.pdf'
        listnamesampdist.append(namesampdist)
        listrtag.append(listpathruns[k])
        listpath.append(pathgmrb)

cmnd = 'mkdir -p ' + pathimag + nameplot
os.system(cmnd)
for path, rtag, namesampdist in zip(listpath, listrtag, listnamesampdist):
    cmnd = 'cp ' + path + ' ' + pathimag + nameplot + '/' + rtag + '_' + namesampdist + '.pdf'
    os.system(cmnd)

