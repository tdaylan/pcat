import os, sys

rtag = sys.argv[1]
pathdata = os.environ["PCAT_DATA_PATH"] + '/data/outp/%s/' % rtag
pathimag = os.environ["PCAT_DATA_PATH"] + '/imag/'
cmnd = 'mkdir -p %s' % pathdata
print cmnd
os.system(cmnd)
cmnd = 'scp -r tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/data/outp/%s/comp.txt %s' % (rtag, pathdata)
print cmnd
os.system(cmnd)
cmnd = 'scp -r tansu@fink2.rc.fas.harvard.edu:/n/fink1/tansu/data/pcat/imag/%s %s' % (rtag, pathimag)
print cmnd
os.system(cmnd)

