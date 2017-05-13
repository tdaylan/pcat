from __init__ import *

listnameplot = [ \
                'post/anim/datacntsene0evt0popA.gif', \
                'post/anim/histdefspop0.gif', \
                'post/anim/histmcutpop0.gif', \
                'post/finl/varbscal/numbpntspop0_trac.pdf', \
                'post/finl/varbscal/fixp_grid0000.pdf', \
                'post/finl/varbscal/fixp_grid0001.pdf', \
                'post/finl/varbscal/fixp_grid0002.pdf', \
                'post/finl/varbscal/fixp_grid0003.pdf', \
                'post/finl/varbscal/fixp_grid0003.pdf', \
                'post/finl/mediconvelemevtApopA.pdf', \
                'post/finl/diag/gmrbmaps_00.pdf', \
                'post/finl/histodim/histdefspop0.pdf', \
                'post/finl/histodim/histdeltllikpop0.pdf', \
                'post/finl/histodim/histdotspop0.pdf', \
                'post/finl/histodim/histdotnpop0.pdf', \
                'post/finl/histodim/histdotmpop0.pdf', \
                'post/finl/histodim/histdotvpop0.pdf', \
                'post/finl/histodim/histmcutpop0.pdf', \
                'post/finl/cmpl/cmpldotspop0.pdf', \
                #'post/finl/cmpl/cmpldefspop0.pdf', \
               ]

listlinegold = [ \
                '20170507_042916_pcat_lens_mock_syst_0002_1000000', \
                '20170507_000608_pcat_lens_mock_syst_0000_1000000', \
                '20170506_213505_pcat_lens_mock_syst_0000_1000000', \
                '20170505_175125_pcat_lens_mock_syst_0004_3000000', \
                '20170505_131822_pcat_lens_mock_syst_0003_1000000', \
                '20170505_095109_pcat_lens_mock_syst_0002_1000000', \
                '20170505_084900_pcat_lens_mock_syst_0002_1000000', \
                '20170504_115457_pcat_lens_mock_perf_0004_1000000', \
                '20170504_071247_pcat_lens_mock_syst_0002_1000000', \
                '20170504_061256_pcat_lens_mock_test_0002_1000000', \
                '20170503_024654_pcat_lens_mock_test_0001_5000000', \
                '20170503_010057_pcat_lens_mock_test_0000_5000000', \
                '20170430_073947_pcat_lens_mock_perf_0005_1000000', \
                '20170430_040715_pcat_lens_mock_perf_0004_1000000', \
                '20170430_000620_pcat_lens_mock_perf_0003_1000000', \
                '20170429_205836_pcat_lens_mock_syst_0002_1000000', \
                '20170429_173210_pcat_lens_mock_syst_0001_1000000', \
                '20170429_122508_pcat_lens_mock_syst_0000_1000000', \
                '20170425_115304_pcat_lens_mock_syst_0002_2000000', \
                '20170429_020452_pcat_lens_mock_perf_0000_5000000', \
                '20170425_083905_pcat_lens_mock_dotn_0002_2000000', \
                '20170425_003054_pcat_lens_mock_syst_0001_2000000', \
                '20170424_223720_pcat_lens_mock_doff_0001_2000000', \
               ]
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
pathlist = pathimag + 'listfink.txt'

#pathgold = pathimag + 'compgold/'
#os.system('rm -rf ' + pathgold)

with open(pathlist) as thisfile:
    listline = thisfile.readlines()
    listline = [x.strip() for x in listline] 

listline = fnmatch.filter(listline, '20*')

print 'compfink initialized...'

for nameplot in listnameplot:   
    
    nameplotshrt = nameplot[:-4].split('/')[-1]
    cmnd = 'mkdir -p ' + pathimag + 'compfink/' + nameplotshrt + '/'
    os.system(cmnd)
    cmnd = 'mkdir -p ' + pathimag + 'compgold/' + nameplotshrt + '/'
    os.system(cmnd)
    print 'Working on ' + cmnd + '...'
    for line in listline:
        
        if line in listlinegold:
            strgtemp = 'gold'
        else:
            strgtemp = 'fink'
        namedest = pathimag + 'comp' + strgtemp + '/' + nameplotshrt + '/' + line + nameplot[-4:]
        if not os.path.isfile(namedest):
            cmnd = 'scp tansu@fink2.rc.fas.harvard.edu:/n/fink2/www/tansu/link/pcat/imag/' + line + '/' + nameplot + ' ' + namedest 
            os.system(cmnd)

