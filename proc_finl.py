import subprocess as subp

def proc_finl():
    
    # post process the samples
    proc_post(gdat)
    
    # make animations
    if gdat.makeanim:
        make_anim(gdat.rtag)

    if gdat.verbtype > 0:
        print 'The ensemble of catalogs is at ' + gdat.pathoutpthis
        if gdat.makeplot:
            print 'The plots are at ' + gdat.pathplot
        print 'PCAT has run successfully. Returning to the OS...'

