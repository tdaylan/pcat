from __init__ import *

for strgextn in ['/imag/', '/data/outp/']:
    
    path = os.environ["PCAT_DATA_PATH"] + strgextn

    for strgfile in os.listdir(path):
        pathfile = path + strgfile
        if os.path.isdir(pathfile) and strgfile[:8].isdigit():
            
            try:
                numbswep = int(pathfile[pathfile.rfind('_')+1:])
            
                #if numbswep < 100000:
                if not os.path.exists(pathfile + '/anim') or not os.listdir(pathfile + '/anim'):
                    os.system('rm -rf ' + pathfile)
            except:
                pass

