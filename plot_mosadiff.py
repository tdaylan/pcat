from __init__ import *
from util import *

rtagfrst = sys.argv[1]
rtagseco = sys.argv[2]

listrtag = [rtagfrst, rtagseco]

indxregiplot = 0

numbcols = 2
indxcols = arange(numbcols)

numbrows = 3
numbsampmosa = numbrows * numbcols

gdatmodi = tdpy.util.gdatstrt()

listgdat = [[] for k in indxcols]
listindxsampmosa = [[] for k in indxcols]
for k in indxcols:
    print 'Processing run tag %s...' % listrtag[k]

    # read gdatfinl objects
    path = os.environ["PCAT_DATA_PATH"] + '/data/outp/' + listrtag[k] + '/gdatfinl'
    try:
        print 'hey'
        listgdat[k] = readfile(path)
    except:
        print 'gdatfinl not found...'
        print
        return

    if numbsampmosa > listgdat[k].numbsamptotl:
        raise Exception('number of samples is less than the number of frames.')
    
    listindxsampmosa[k] = choice(listgdat[k].indxsamptotl, size=numbsampmosa, replace=False)

# common gdat object
gdat = listgdat[0]

for d in gdat.indxregi:
    for i in gdat.indxener:
        for m in gdat.indxevttplot:
            figr, axgr = plt.subplots(numbrows, numbcols, figsize=(numbcols * gdat.plotsize, numbrows * gdat.plotsize))
            for a, axrw in enumerate(axgr):
                for b, axis in enumerate(axrw):
                    
                    n = listindxsampmosa[b][a]
                    gdatmodi.thissampvarb = listgdat[b].listsampvarb[n, :].flatten()
                    gdatmodi.thissamp = listgdat[b].listsamp[n, :].flatten()
                    
                    if listgdat[b].fittnumbtrap > 0:
                        gdatmodi.thisindxelemfull = listgdat[b].listindxelemfull[n]
                        proc_samp(listgdat[b], gdatmodi, 'this', 'fitt')

                    if a == numbrows - 1:
                        axis.set_xlabel(gdat.labllgaltotl)
                    else:
                        axis.set_xticklabels([])
                    if b == 0:
                        axis.set_ylabel(gdat.lablbgaltotl)
                    else:
                        axis.set_yticklabels([])
                        
                    imag = retr_imag(gdat, axis, gdat.cntpdata[indxregiplot], '', 'fitt', 'cntpdata', i, m)
                    for l in listgdat[k].fittindxpopl:
                        supr_fram(gdat, gdatmodi, 'this', 'fitt', axis, indxregiplot, l)
                
            if gdat.enerbins:
                plt.figtext(0.5, 0.93, gdat.strgener[i], ha='center', va='center')
            axiscomm = figr.add_axes([0.92, 0.1, 0.02, 0.8])
            cbar = figr.colorbar(imag, cax=axiscomm)
            cbar.set_ticks(gdat.tickcntpdata)
            cbar.set_ticklabels(gdat.lablcntpdata)
            plt.subplots_adjust(left=0.1, top=.91, hspace=0.03, wspace=0.1, bottom=0.09)
            if l == 1:
                strg = ''
            else:
                strg = 'pop%d' % l
            pathfinl = getattr(gdat, 'path' + gdat.namesampdist + 'finl')
            if m == None:
                path = pathfinl + 'mosa' + strg + 'reg%dene%dA.pdf' % (indxregiplot, gdat.indxenerincl[i])
            else:
                path = pathfinl + 'mosa' + strg + 'reg%dene%devtt%d.pdf' % (indxregiplot, gdat.indxenerincl[i], gdat.indxevttincl[m])
            make_legdmaps(gdat, 'this', 'fitt', axis)
            print 'Writing to %s...' % path
            figr.savefig(path)
            plt.close(figr)


