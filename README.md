# PCAT (Probabilistic Cataloger)

PCAT is a Bayesian framework to sample from the catalog space. 

## Usage

### Options
`verbtype`
Verbosity level
- 0: no standard output
- 1: print only critical statements along with periodic progress
- 2: full log of internal variables, to be used only for diagnostic purposes
Default: 1

`plotperd`
Number of sweeps between frames
Default: 50000

`makeplot`
Boolean flag to allow making of plots
Default: True

`di`


         diagsamp=False, \
         numbswep=2000000, \
         numbburn=None, \
         factthin=None, \
         regitype='ngal', \
         datatype='inpt', \
         randinit=True, \
         maxmgang=None, \
         minmspec=None, \
         maxmspec=None, \
         minmsind=None, \
         maxmsind=None, \
         meansind=None, \
         stdvsind=None, \
         minmfdfnnorm=None, \
         maxmfdfnnorm=None, \
         minmfdfnslop=None, \
         maxmfdfnslop=None, \
         fdfntype='powr', \
         psfntype=None, \
         proppsfn=True, \
         numbpopl=1, \
         indxevttincl=arange(2, 4), \
         indxenerincl=arange(5), \
         maxmnumbpnts=array([1000]), \
         initnumbpnts=None, \
         trueinfo=False, \
         pntscntr=False, \
         numbproc=None, \
         liketype='pois', \
         pixltype='heal', \
         exprtype='ferm', \
         lgalcntr=0., \
         bgalcntr=0., \
         margsize=None, \
         maxmnormback=None, \
         minmnormback=None, \
         numbsidecart=None, \
         numbsideheal=None, \
         maxmangleval=None, \
         spmrlbhl=2., \
         stdvfdfnnorm=0.05, \
         stdvfdfnslop=0.1, \
         stdvpsfipara=0.1, \
         stdvback=0.04, \
         stdvlbhl=0.1, \
         stdvspec=0.15, \
         stdvpropsind=0.15, \
         fracrand=0.05, \
         mocknumbpnts=None, \
         mockfdfnslop=None, \
         mockfdfnsloplowr=None, \
         mockfdfnslopuppr=None, \
         mockfdfnbrek=None, \
         mocknormback=None, \
         mockfdfntype='powr', \
         mockpsfntype=None, \
         strgexpr=None, \
         strgback=None, \
         lablback=None, \
         nameback=None, \
         strgexpo=None, \
         probprop=None, \
        ):
    
    # start the timer
    timetotlreal = time.time()
    timetotlproc = time.clock()
   
    # defaults
    ## Fermi-LAT
    if exprtype == 'ferm':
        if maxmgang == None:
            maxmgang = 20.
        if minmsind == None:
            minmsind = array([1.2])
        if maxmsind == None:
            maxmsind = array([3.2])
        if meansind == None:
            meansind = array([2.2])
        if stdvsind == None:
            stdvsind = array([0.3])
        if minmfdfnnorm == None:
            minmfdfnnorm = array([1e0])
        if maxmfdfnnorm == None:
            maxmfdfnnorm = array([1e2])

