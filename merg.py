from __init__ import *
from PyPDF2 import PdfFileMerger

pdfs = ['file1.pdf', 'file2.pdf', 'file3.pdf', 'file4.pdf']

pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
listpathruns = fnmatch.filter(os.listdir(pathimag), '20*')
numbruns = len(listpathruns)
listpath = []
for k in range(numbruns):
    for namesampdist in ['prio', 'post']:
        pathgmrb = pathimag + listpathruns[k] + '/' + namesampdist + '/finl/diag/gmrbmaps_00.pdf'
        listpath.append(pathgmrb)

merger = PdfFileMerger()

for path in listpath:
    if os.path.exists(path):
        print path
        merger.append(path)

#merger.write("result.pdf")

