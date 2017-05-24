from __init__ import *

#listnameplot = [ \
#                'post/anim/datacntsene0evt0popA.gif', \
#                'post/anim/histdefspop0.gif', \
#                'post/anim/histmcutpop0.gif', \
#                'post/finl/varbscal/lgalsour_trac.pdf', \
#                'post/finl/varbscal/numbpntspop0_trac.pdf', \
#                'post/finl/varbscal/fracsubh_trac.pdf', \
#                'post/finl/varbscal/fixp_grid0000.pdf', \
#                'post/finl/varbscal/fixp_grid0001.pdf', \
#                'post/finl/varbscal/fixp_grid0002.pdf', \
#                'post/finl/varbscal/fixp_grid0003.pdf', \
#                'post/finl/varbscal/fixp_grid0003.pdf', \
#                'post/finl/mediconvelemevtApopA.pdf', \
#                'post/finl/diag/gmrbmaps_00.pdf', \
#                'post/finl/histodim/histdefspop0.pdf', \
#                'post/finl/histodim/histdeltllikpop0.pdf', \
#                'post/finl/histodim/histrelepop0.pdf', \
#                'post/finl/histodim/histrelnpop0.pdf', \
#                'post/finl/histodim/histreldpop0.pdf', \
#                'post/finl/histodim/histrelcpop0.pdf', \
#                'post/finl/histodim/histmcutpop0.pdf', \
#                'post/finl/llik_hist.pdf', \
#                'post/finl/llik_trac.pdf', \
#                'post/finl/histdefl.pdf', \
#                'post/finl/histdeflelem.pdf', \
#                'post/finl/mosapop0ene0evtt0.pdf', \
#                'post/finl/medideflcompevtApopA.pdf', \
#                'post/finl/medimagnresievtApopA.pdf', \
#                'post/finl/medimagnpercresievtApopA.pdf', \
#                'post/finl/postdeflsing0000popA.pdf', \
#                'post/finl/postdeflsing0001popA.pdf', \
#                'post/finl/postdeflsing0002popA.pdf', \
#                'post/finl/postdeflsing0003popA.pdf', \
#                'post/finl/postdeflsingresi0000popA.pdf', \
#                'post/finl/postdeflsingresi0001popA.pdf', \
#                'post/finl/postdeflsingresi0002popA.pdf', \
#                'post/finl/postdeflsingresi0003popA.pdf', \
#                'post/finl/cmpl/cmpldotspop0.pdf', \
#                #'post/finl/cmpl/cmpldefspop0.pdf', \
#               ]
#
listlinegold = [ \
                '20170522_020917_pcat_lens_mock_syst_0001_5000000', \
                '20170521_232045_pcat_lens_mock_syst_0000_5000000', \
               ]
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'

#cmnd = 'ssh tansu@fink2.rc.fas.harvard:ls /n/fink1/tansu/data/pcat/imag/ | xargs -n 1 basename > /n/fink1/tansu/data/pcat/imag/listfink.txt'
#os.system(cmnd)
cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/imag/listfink.txt %s' % pathimag
os.system(cmnd)
pathlist = pathimag + 'listfink.txt'

pathgold = pathimag + 'compgold/'
os.system('rm -rf ' + pathgold)

with open(pathlist) as thisfile:
    listline = thisfile.readlines()
    listline = [x.strip() for x in listline] 

#strgsrch = '20*'
strgsrch = '2017052*'
listline = fnmatch.filter(listline, strgsrch)

print 'compfink initialized...'

def comp(nameplot):
    
    if nameplot == 'DS_Store':
        return

    nameplotshrt = nameplot[:-4].split('/')[-1]
    cmnd = 'mkdir -p ' + pathimag + 'compfink/' + nameplotshrt + '/'
    #os.system(cmnd)
    print cmnd
    cmnd = 'mkdir -p ' + pathimag + 'compgold/' + nameplotshrt + '/'
    #os.system(cmnd)
    print cmnd
    print 'Working on ' + nameplotshrt + '...'
    for line in listline:
        
        if line in listlinegold:
            strgtemp = 'gold'
        else:
            strgtemp = 'fink'
        namedest = pathimag + 'comp' + strgtemp + '/' + nameplotshrt + '/' + line + nameplot[-4:]
        if not os.path.isfile(namedest):
            cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink2/www/tansu/link/pcat/imag/' + line + '/' + nameplot + ' ' + namedest 
            #os.system(cmnd)
            print cmnd
            print
    print
    print
    print
    print
    print

pathrefr = os.environ["PCAT_DATA_PATH"] + '/imag/20170521_232130_pcat_lens_mock_syst_0000_2000000/'
for namefrst in os.listdir(pathrefr):
    namefrstprim = os.path.join(pathrefr, namefrst)
    if os.path.isdir(namefrstprim):
        for nameseco in os.listdir(namefrstprim):
            namesecoprim = os.path.join(pathrefr, nameseco)
            if os.path.isdir(namesecoprim):
                for namethrd in os.listdir(namesecoprim):
                    namethrdprim = os.path.join(pathrefr, namethrd)
                    comp(namethrdprim)
            elif nameseco.endswith('pdf') or nameseco.endswith('gif'):
                comp(namesecoprim)
    elif namefrst.endswith('pdf') or namefrst.endswith('gif'):
        comp(namefrstprim)



