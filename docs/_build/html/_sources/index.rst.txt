Welcome to PCAT's documentation!
================================

When testing hypotheses or inferring their free parameters, a recurring problem is to compare models that contain a number of objects whose multiplicity is itself unknown. Therefore, given some data, it is desirable to be able to compare models with different number of parameters. One way of achieving this is to obtain a point estimate (usually the most likely point) in the parameter space of each model and, then, rely on some information criterion to penalize more complex models for excess degrees of freedom. Another way is to sample from the parameter space of each model and compare their Bayesian evidence. Yet another is to take samples from the union of the hypotheses using a set of transdimensional jumps across models. This is what PCAT (Probabilistic Cataloger) is designed for.

PCAT is a hierarchical, transdimensional MCMC sampler. Given some data, it can be used to sample from the **catalog space**, i.e., the hypothesis space of a model with a variable number of components. Alternatively, it can be used as a general purpose Poisson mixture sampler to infer clusters in a dataset.

A user manual for PCAT will be published here soon.

.. _sectinst:

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

.. function:: pcat.main.init(...)

    Given an observed dataset, sample from the hypothesis space.

    **Sampler settings**

    :param numbswep: Number of samples to be taken by each process

    :type numbswep: int


    :param numbburn: Number of samples to be discarded from the beginning of each chain

    :type numbburn: int


    :param factthin: Factor by which to thin each chain. Only one sample out of ``factthin`` samples is saved.

    :type factthin: int


    :param numbproc: Number of processes. The total number of samples before thinning and burn-in, will be ``numbproc`` times ``numbswep``.

    :type numbproc: int


    **Input**
    :param indxenerincl: Indices of energy bins to be taken into account. It is only effective if data is provided by the user, i.e., for non-simulation runs, where it defaults to all available energy bins. Other energy bins are discarded.
    
    :type indxenerincl: ndarray int

    
    :param indxevttincl: Indices of PSF class bins to be taken into account. Works similar to ``indxenerincl``.
    
    :type indxevttincl: ndarray int

   
    **Output**

    :param verbtype: Verbosity level
        
        - ``0`` No standard output
        - ``1`` Minimal standard output including status of progress (Default)
        - ``2`` Diagnostic verbose standard output
    
    :type verbtype: int

    
    :param pathbase: Data path of PCAT. See :ref:`sectoutp`.

    :type pathbase: str


    **Asscociations with the reference catalog**
 
    :param anglassc: Radius of the circle within which sample catalog point sources can be associated with the point sources in the reference catalog.
    
    :type anglassc: float


    :param margfactcomp: The ratio of the side of the square in which the sample point sources are associated with the reference catalog, to the size of the image.

    :type margfactcomp: float


    :param nameexpr: A string that describes the provided reference catalog to be shown in the plot legends.

    :type nameexpr: str
    

    :param cntrpnts: Force the mock data to a single PS at the center of the image. Defaults to ``False``.
    
    :type cntrpnts: bool


    **Unfunctional**
 
    :param numbspatdims: Number of spatial dimensions. Currently not functional.

    :type numbspatdims: None


    **Asscociations with the reference catalog**

    :param pntstype: Functional type of point sources. 
        
        - ``lght`` Point sources are light sources
        - ``lens`` Point sources are lenses.
    
    :type pntstype: str


    :param randinit: Force the initial state to be randomly drawn from the prior. Default behavior for mock data is to initialize the chain with the true state.

    :type randinit: bool


    **Adaptive burn-in**
    :param loadvaripara: Load the diagonal elements of the previously learned covariance

    :type loadvaripara: None


    :param optiprop: Optimize the scale of each proposal by acceptance rate feedback. All samples during the tuning are discarded.
        
    :type randinit: bool


    **Post processing**
    :param regulevi: Regularize the log-evidence estimate.

    :type regulevi: bool


    :param strgexprflux: Name of the FITS file (without the extension) in ``pathdata`` containing the observed data as an ``ndarray``. The file should contain a numpy array of dimension $(N_e, N_{pix}, N_{psf}$, where $N_e$ is the number of energy bins, $N_{pix}$ is the number of spatial pixels and $N_{psf}$ is the number of PSF classes. The units should be photons per $cm^2$ per seconds per GeV.

    :type strgexprflux: str


    :param strgcatl: A descriptive name for the provided reference catalog to be shown in the plot legends.

    :type strgcatl: str


    :param strgback: A list of FITS file names (without the extension) in ``pathdata`` each containing a spatial template for the background prediction as an ``ndarray``. See ``strgexprflux`` for the content of the file and its unit. One element of the list can be a float, indicating an isotropic template with the provided amplitude.

    :type strgback: list of str or int
    

    :param lablback: a list of axis labels for the spatial background templates to be shown in plots.

    :type lablback: list of str


    :param strgexpo: Name of the FITS file (without the extension) in ``pathdata`` containing the exposure map. See ``strgexprflux`` for the format of the numpy array. ``strgexpo`` can also be a float, in which case the exposure map will be assumed to be uniform across along all data dimensions.

    :type strgexpo: str or float


    :param liketype: Type of the likelihood. 

        - ``pois`` Poisson probability of getting the observed number of counts given the model prediction (default). 
        - ``gaus`` Gaussian approximation of the above. This may accelerate the execution in cases, where the bottle neck of the sampler time budget is likelihood evaluation.
    
    :type liketype: strg


    :param exprtype: Name of the experiment used to collect the observed data. ``exprtype`` can be used to set other options to their default values for the particular experiment. 
        - ``ferm`` Fermi-LAT
        - ``chan`` Chandra
        - ``hubb`` HST
        - ``sdss`` SDSS

    :type exprtype: str


    :param lgalcntr: Galactic longitude of the image center. ``lgalcntr`` and ``bgalcntr`` are used to rotate the observed data, exposure and background maps as well as the provided reference catalog to the center of the ROI. They are only effective when pixelization is HealPix, i.e, ``pixltype='heal'``.

    :type lgalcntr: float


    :param bgalcntr: Galactic latitude of the image center. See ``lgalcntr``.

    :type bgalcntr: float

 
    :param maxmangl: Maximum angular separation at which PSF can be interpolated. It defaults to three times the diagonal legth of the image, enough to evaluate the PSF across the whole image.
   
    :type maxmangl: float


    :param pixltype: Type of the pixelization.
        
        - ``heal`` HealPix
        - ``chan`` Cartesian
        - ``unbd`` Unbinned (Not yet functional)


    **Plotting**

    :param makeplot: Make output plots, which is the default behavior. If ``False``, no output plots are produced.

    :type makeplot: bool


    :param numbswepplot: Number of samples (before thinning and burn-in) for which one set of frame plots will be produced. Frame plots reveal individual samples in detail and are later used for producing animations.

    :type numbswepplot: int


    :param scalmaps: A string that sets the stretch of the count maps

        - ``asnh`` Arcsinh (default)
        - ``self`` Linear
        - ``logt`` Log 10

    :type scalmaps: str
    

    :param satumaps: Saturate the count maps

    :type satumaps: bool


    :param exprinfo: Overplot the provided reference catalog on the output plots.

    :type exprinfo: bool


    :param makeanim: Make animations of the frame plots. Defaults to ``True``.

    :type makeanim: bool


    :param anotcatl: Anotate the catalog members on the plots, if an annotation text is provided along with the reference catalog. (Default: ``False``)

    :type anotcatl: bool

    
    :param strgbinsener: A string holding the label for the energy axis.

    :type strgbinsener: str


    **Diagnostics**

    :param diagmode: Start the run in diagnostic mode. Defaults to ``False``.

    :type diagmode: bool


    :param strgexprname: A string describing the experiment used to collect the observed data.

    :type strgexprname: str


    :param strganglunit: Label for the spatial axes.

    :type strganglunit: str


..    
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

..
.. Features
.. ----------
.. 
.. Transdimensionality
.. +++++++++++++++++++
.. 
.. Hierarchical priors
.. +++++++++++++++++++++
.. 
.. Adaptive burn-in
.. +++++++++++++++++
.. 
.. Labeling degeneracy
.. ++++++++++++++++++++++
.. 
.. Input
.. --------------
.. 
.. Producing mock data
.. ++++++++++++++++++++++
.. 
.. Prior and initial state specification
.. ++++++++++++++++++++++
.. 

.. _sectoutp:

Output
-------------

PCAT expects input data in ``pathbase/data/inpt/`` and writes the output chains and plots in the relevant subfolders ``pathbase/data/outp/rtag`` and ``pathbase/imag/rtag``, respectively, where ``rtag`` is the run tag. The ``pathbase`` folder is created if it does not already exist. It defaults to the value of the environment variable ``PCAT_DATA_PATH``. Therefore ``pathbase`` can be omitted by setting the environment variable ``PCAT_DATA_PATH``.

Plots
+++++
If not disabled by the user, PCAT produces plots in every stage of a run. Some plots are produced in the initial setup, frame plots are produced at predetermined times during the sampling and others are produced in the postprocessing after all chains have run. There are five subfolders in the plot folder of a given run.

- ``init`` Initial setup, pixel lookup table (if applicable)
- ``fram`` Frame plots
- ``diag`` Diagnostic plots
- ``post`` Posterior distribution plots of model parameters and derived quantities, prior and likelihood.
- ``anim`` GIF animations made from the frame plots in ``fram`` that are producing during sampling.

.. Chain
.. +++++


.. Tutorial
.. --------------
.. In this tutorial, a typical run on a mock dataset will be illustrated. Assuming that you have :ref:`installed <sectinst>`

.. Diagnosing the sampler
.. --------------
.. .note:: The acceptance rate of the birth and death moves shows whether the prior on the amplitude of the elements is appropriate. If the acceptance rate of birth and deaths proposals is too low, the lower limit of the prior on the amplitude is too high. Likewise, if it is too low, the lower limit of the prior on the amplitude is too low. This behaviour is due to the fact that element parameters (position, amplitude and, if relevant, spectral or shape parameters) are drawn randomly from the prior when a birth is proposed. Therefore the lower limit on the amplitude prior should be adjusted such that the birth and death acceptance rate is between $\sim 5%% - \sim30%%$. Otherwise across-model sampling will be inefficient and result in slow convergence.  
..
.. Post processing samples from the catalog space
.. ++++++++++++++++++++
.. 
.. 


Diagnostics
------------

.. note::

   Make sure to run PCAT with ``diagmode=False`` for good time performance. ``diagmode=True`` puts the sampler in a conservative diagnostic mode and performs extensive checks on the state of critical data structures to ensure that the model and proposals are self consistent.

.. toctree::
    :maxdepth: 2


Garbage collection
-------------------

PCAT produces two folders for the output of each run, one for plots and the other to contain the chain saved to the disc. Given that many test and intermediate runs may be needed before each science run, the number of folders (and files therein) may increase quickly. In order to avoid this, the script ``gcol.py`` is provided to the user as a convenience. When executed, it erases the output from all runs that

- have not run to completion, i.e., does not have the animations, which are produced at the very end, or
- has collected less than 100000 samples per chain.

