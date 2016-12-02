Welcome to PCAT's documentation!
================================

When testing hypotheses or inferring their free parameters, a recurring problem is to compare models that contain a number of objects whose multiplicity is itself unknown. Therefore, given some data, it is desirable to be able to compare models with different number of parameters. One way of achieving this is to obtain a point estimate (usually the most likely point) in the parameter space of each model and, then, rely on some information criterion to penalize more complex models for excess degrees of freedom. Another way is to sample from the parameter space of each model and compare their Bayesian evidence. Yet another is to take samples from the union of the hypotheses using a set of transdimensional jumps across models. This is what PCAT (Probabilistic Cataloger) is designed for.

PCAT is a hierarchical, transdimensional MCMC sampler. Given some data, it can be used to sample from the **catalog space**, i.e., the hypothesis space of a model with a variable number of components. Alternatively, it can be used as a general purpose Poisson mixture sampler to infer clusters in a dataset.

A user manual for PCAT will be published here soon.

Installation
-------------
To install PCAT you can use pip

.. code::

    pip install pcat

or download `the latest release <https://github.com/tdaylan/pcat/releases/>`_ and run

.. code::

    python setup.py install


Usage
---------

All user interaction with PCAT is accomplished through the ``pcat.main.init()`` function.

.. function:: pcat.main.init()

    Given an observed dataset, sample from the hypothesis space.

    # User Interaction
    
    :param verbtype:

        Verbosity level
        
        - ``0`` No standard output
        - ``1`` Minimal standard output including status of progress (Default)
        - ``2`` Diagnostic verbose standard output
    
    :type verbtype: int

    
    :param pathbase: Data path of PCAT. PCAT expects input data in ``pathbase/data/inpt/`` and writes the output chains and plots in the relevant subfolders ``pathbase/data/outp/rtag`` and ``pathbase/imag/rtag``, respectively, where ``rtag`` is the run tag. The ``pathbase`` folder is created if it does not already exist. It defaults to the value of the environment variable ``PCAT_DATA_PATH``. Therefore ``pathbase`` can be omitted by setting the environment variable ``PCAT_DATA_PATH``.

    :type pathbase: str



    # diagnostics

    :param diagmode:
        
        Start the run in diagnostic mode. Defaults to ``False``.

    :type diagmode: bool


    :param cntrpnts:
        
        Force the mock data to a single PS at the center of the image. Defaults to ``False``.
    
    :type cntrpnts: bool


    # chain
    
    :param numbswep:
        
        Number of samples to be taken by each process

    :type numbswep: int


    :param numbburn:
        
        Number of samples to be discarded from the beginning of each chain

    :type numbburn: int


    :param factthin:
        
        Factor by which to thin each chain. Only one sample out of ``factthin`` samples is saved.

    :type factthin: int


    :param indxenerincl:
        
        Indices of energy bins to be taken into account. It is only effective if data is provided by the user, i.e., for non-simulation runs, where it defaults to all available energy bins. Other energy bins are discarded.
    
    :type indxenerincl: ndarray int

    
    :param indxevttincl:
        
        Indices of PSF class bins to be taken into account. Works similar to ``indxenerincl``.
    
    :type indxevttincl: ndarray int

    
         
..     
         # comparison with the reference catalog
         anglassc=None, \
         margfactcomp=0.9, \
         nameexpr=None, \
         numbspatdims=2, \
         pntstype='lght', \
         randinit=None, \
         loadvaripara=False, \
         optiprop=False, \
         regulevi=False, \
         strgexpr=None, \
         strgcatl=None, \
         strgback=[1.], \
         lablback=None, \
         nameback=None, \
         strgexpo=1., \
         numbproc=None, \
         liketype='pois', \
         exprtype='ferm', \
         lgalcntr=0., \
         bgalcntr=0., \
         maxmangl=None, \
         exprinfo=None, \
         pixltype=None, \
         
         # plotting
         numbswepplot=50000, \
         makeplot=True, \
         scalmaps='asnh', \
         satumaps=None, \
         makeanim=True, \
         anotcatl=False, \
         strgbinsener=None, \
         strgexprname=None, \
         strganglunit=None, \
         strganglunittext=None, \
         anglfact=None, \
         fluxfactplot=None, \
         enerfact=None, \
         
         # misc
         strgfunctime='clck', \
         strgxaxi=None, \
         strgyaxi=None, \
         # model
         ## PSF
         specfraceval=0.1, \
         numbangl=1000, \
         binsangltype='logt', \
         numbsidepntsprob=400, \
         strgfluxunit=None, \
         strgflux=None, \
         strgenerunit=None, \
         indxenerfull=None, \
         indxevttfull=None, \
         binsenerfull=None, \
         maxmnumbpnts=array([1000]), \
         asymfluxprop=False, \
         psfninfoprio=True, \
         ## spectral
         # prior
         priotype='logt', \
         priofactdoff=0., \
         margfactmodl=0.9, \
         bindprio=False, \
         maxmbacp=None, \
         minmbacp=None, \
         maxmgang=None, \
         minmmeanpnts=None, \
         maxmmeanpnts=None, \
    
         spatdisttype=None, \
         spatdistslop=None, \
         
         fluxdisttype=None, \
         minmfluxdistslop=None, \
         maxmfluxdistslop=None, \
         minmfluxbrek=None, \
         maxmfluxbrek=None, \
         minmfluxdistbrek=None, \
         maxmfluxdistbrek=None, \
         minmfluxdistsloplowr=None, \
         maxmfluxdistsloplowr=None, \
         minmfluxdistslopuppr=None, \
         maxmfluxdistslopuppr=None, \
         minmsinddistmean=None, \
         maxmsinddistmean=None, \
         minmsinddiststdv=None, \
         maxmsinddiststdv=None, \
         psfntype=None, \
         varioaxi=None, \
         minmsigm=None, \
         maxmsigm=None, \
         meansigm=None, \
         stdvsigm=None, \
         minmgamm=None, \
         maxmgamm=None, \
         meangamm=None, \
         stdvgamm=None, \
         minmpsff=None, \
         maxmpsff=None, \
         meanpsff=None, \
         stdvpsff=None, \
         minmfluxsour=None, \
         maxmfluxsour=None, \
         minmsizesour=None, \
         maxmsizesour=None, \
         minmratisour=None, \
         maxmratisour=None, \
         minmellphost=None, \
         maxmellphost=None, \
         minmsherhost=None, \
         maxmsherhost=None, \
         minmbeinhost=None, \
         maxmbeinhost=None, \
    
         spectype=None, \
         
         curvdistmean=None, \
         curvdiststdv=None, \
         
         minmflux=None, \
         maxmflux=None, \
        
         # proposals
         numbpntsmodi=1, \
         stdvprophypr=0.1, \
         stdvproppsfp=0.01, \
         stdvpropbacp=0.01, \
         stdvproplenp=0.01, \
         stdvlbhl=0.1, \
         stdvlbhlvari=True, \
         stdvflux=0.15, \
         stdvspep=0.15, \
         stdvspmrsind=0.2, \
         probrand=0.05, \
         boolpropfluxdist=True, \
         boolpropfluxdistbrek=True, \
         prophypr=True, \
         proppsfp=True, \
         propbacp=True, \
         proplenp=True, \
         propcomp=True, \
         probtran=None, \
         probbrde=1., \
         radispmr=None, \
         truevarioaxi=None, \
         truepsfntype=None, \
         # mock data
         mockspatdisttype=None, \
         mockspatdistslop=None, \
         mockfluxdisttype=None, \
         mockminmflux=None, \
         mockmaxmflux=None, \
         mockfluxdistslop=None, \
         mockfluxdistbrek=None, \
         mockfluxdistsloplowr=None, \
         mockfluxdistslopuppr=None, \
         mockspectype=None, \
         mocksinddistmean=None, \
         mocksinddiststdv=None, \
         mockpsfntype=None, \
         mockvarioaxi=None, \
         mockstrgback=None, \
         mockbacp=None, \
         
         mocklgalsour=None, \
         mockbgalsour=None, \
         mockfluxsour=None, \
         mocksizesour=None, \
         mockratisour=None, \
         mockanglsour=None, \
         mocklgalhost=None, \
         mockbgalhost=None, \
         mockellphost=None, \
         mockanglhost=None, \
         mocksherhost=None, \
         mocksanghost=None, \
         mockbeinhost=None, \
         
         mocknumbpnts=None, \
         numbsidecart=200, \
         numbsideheal=256, \
         numbdatasamp=100, \


Features
----------

Transdimensionality
+++++++++++++++++++

Hierarchical priors
+++++++++++++++++++++

Adaptive burn-in
+++++++++++++++++

Initial state
++++++++++++++++++++++

Identity of components
++++++++++++++++++++++

.. toctree::
    :maxdepth: 2


