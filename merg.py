from __init__ import *
from PyPDF2 import PdfFileMerger

pdfs = ['file1.pdf', 'file2.pdf', 'file3.pdf', 'file4.pdf']

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

merger = PdfFileMerger()

cmnd = 'mkdir -p ' + pathimag + nameplot
print cmnd
#os.system(cmnd)
for path, rtag, namesampdist in zip(listpath, listrtag, listnamesampdist):
    cmnd = 'cp ' + path + ' ' + pathimag + nameplot + '/' + rtag + '_' + namesampdist + '.pdf'
    print cmnd
    os.system(cmnd)
    #if os.path.exists(path):
    #    merger.append(path)
#merger.write(pathimag + nameplot)

