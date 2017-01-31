Welcome to PCAT's documentation!
================================

When testing hypotheses or inferring their free parameters, a recurring problem is to compare models that contain a number of objects whose multiplicity is itself unknown. Therefore, given some data, it is desirable to be able to compare models with different number of parameters. One way of achieving this is to obtain a point estimate (usually the most likely point) in the parameter space of each model and, then, rely on some information criterion to penalize more complex models for excess degrees of freedom. Another way is to sample from the parameter space of each model and compare their Bayesian evidence. Yet another is to take samples from the union of the hypotheses using a set of transdimensional jumps across models. This is what PCAT (Probabilistic Cataloger) is designed for.

PCAT is a hierarchical, transdimensional MCMC sampler. It's theoretical framework is introduced in `Daylan, Portillo & Finkbeiner (2016) <https://arxiv.org/abs/1607.04637>`_, submitted to ApJ. Given some data, it can be used to sample from the **catalog space**, i.e., the hypothesis space of a model with a variable number of components. Alternatively, it can be used as a general purpose Poisson mixture sampler to infer clusters in a dataset.

.. toctree::
   :maxdepth: 4


.. _sectinst:


Installation
-------------
To install PCAT you can use pip

.. code::

    pip install pcat

or download `the latest release <https://github.com/tdaylan/pcat/releases/>`_ and run

.. code::

    python setup.py install


.. _sectinpt:


Features
----------
Compared to mainstream point source inference methods, PCAT has a series of desirable features.

- samples from the space of catalogs given some observation, unlike conventional cataloging, which estimates the most likely catalog,
- allows marginalization over all relevant nuisance parameters in the problem, including the dimensionality of the nuisance.
- reveals potentially non-Gaussian within and across model covariances
- constrains source population characteristics via hierarchical priors,
- is a Bayesian framework, because point estimates fail in nearly degenerate likelihood topologies
- implements Occam's razor, i.e., model parsimony, through detailed balance across models, 
- does not discard information contained in low-significance ($< 4 \sigma$) fluctuations in the observed dataset,
- reduces to a deterministic cataloger when the labeling degeneracy is explicitly broken,
- simultaneously infers the PSF and the level of diffuse background.

Transdimensionality
+++++++++++++++++++

PCAT takes steps across models by adding parameters drawn from the prior or killing them while maintaining detailed balance in the hyper model space.


Hierarchical priors
+++++++++++++++++++++

When there are multiple model sources, each with a set of properties, it is more natural to put priors on the distribution of these properties, as opposed to placing individual priors separately on each source property. This assumes that a certain property of all point sources in a given population are drawn from a single probability distribution. This is particularly useful when such a parametrization is subject to inference, where individual sources can be marginalized over. This results in a hierarchical prior structure, where the prior is placed on the distribution of source properties, i.e., **hyperparameters**, and the prior on the individual source properties are made conditional on these hyperparameters. 


.. Proposals
.. +++++++++++++++++++++


Adaptive burn-in
+++++++++++++++++

PCAT discards the first ``numbburn`` samples and thins the resulting chain by a factor ``factthin``. In order to estimate proposal scales for parameters, it initially takes the derivative of log-likelihood with respect to all parameters at a likely region in the parameter space.

Labeling degeneracy
++++++++++++++++++++++

Due to **hairlessness** of the point sources, the likelihood function is invariant to their permutations in the parameter vector, i.e., exchanging the labels of two point sources leaves the likelihood invariant. This fact has consequences for **nonpersistent** sources, which get born or killed at least once during an MCMC run. Because of label changes, the posterior of these parameters look the same, which makes them useless for inferring their properties. In order to constrain such sources, the degeneracy must be broken in post-processing of the samples. Note that, if the sampling is continued sufficiently long, e.g., for a Hubble time, the posteior of all transdimensional parameters will eventually look similar.

.. Breaking the labeling degeneracy
.. +++++++++++++++++++

Performance
+++++++++++++++++++

The above properties are made possible by enlarging the hypothesis space so much that there is no mismodeling of the observed data. This is, however, at the expense of seriously slowing down inference. The main reason for thfthese are possible at the expense of PCAT being a much slower algorithm.


Input
--------------
PCAT expects input data in the folder ``pathbase/data/inpt/``. The principal input is the observed dataset, which could be binned or unbinned (currently unfunctional).

Supplying data
+++++++++++++++
Input dataset is provided through the ``strgexprflux`` argument. This should be the name of a FITS file (including the ``.fits`` extension), which contains a numpy array of dimension :math:`N_e \times N_{pix} \times N_{psf}`, where :math:`N_e` is the number of energy bins, :math:`N_{pix}` is the number of spatial pixels and :math:`N_{psf}` is the number of PSF classes. The units should be photons per cm :math:`^2` per seconds per GeV.

The exposure map is suppled via the ``strgexpo`` argument. This should be the name of a FITS file (including the ``.fits`` extension). The format should be the same as ``strgexprflux``, whereas units should be cm :math:`^2` GeV.

Similary, background templates can be provided via the ``back`` argument. ``back`` should be a list of FITS file names (including the ``.fits`` extension), whose format should be same as that of ``strgexprflux``.

Specifying model priors
+++++++++++++++++++++++++

The prior structure of the model is set by the relevant arguments to ``pcat.main.init()``. PCAT allows the following components in the model:

- Background prediction

    Given some an observed dataset, it can be expressed as a Poisson realization of the integral emission from a large number of point sources. However, due to the prohibitively large number of point sources required, it is favorable to represent the contribution of extremely faint point sources (those that negligibly affect the likelihood) with background templates. PCAT allows a list of spatially and spectrally distinct background templates to be in the model simultaneously. 


- Point Spread Function (PSF)
    
    PSF defines a how a delta function in the position space projects onto the data space, e.g., the image. PCAT assumes that all point sources are characterized by the same PSF. 


Generating mock data
+++++++++++++++++++++
PCAT ships with a built-in mock (simulated) data generator. The generated mock data is randomly drawn from a generative mock model, not to be confused with the above model subject to inference. Once the user configures the prior probability density of this model, PCAT samples from the catalog space given the simulated dataset. Some of the mock model parameters can be fixed, i.e., assigned delta function priors. These are

- the number of mock point sources, ``mocknumbpnts``,
- hyperparameters controlling the population characteristics of these point sources.

All other mock model parameters are fair draws from the hierarchical prior. 

When working on variations of a certain problem or different analyses on the same dataset, it is useful to have default priors. PCAT allows unique defaults for different built-in experiments, controlled by the argument ``exprtype``. Currently the built-in experimental types are 

- ``ferm``: Fermi-LAT (Default)
- ``chan``: Chandra
- ``hubb``: Hubble Space Telescope
- ``sdss``: SDSS
- ``gaia``: Gaia

By setting ``exprtype`` the user imposes the default prior structure for the chosen experimental type. However, it is also desirable to be able to change the prior of a specific parameter. This is accomplished by setting the :ref:`relevant argument(s) <sectoutp>`.

.. note::

    The mock generative model defaults to that of the model subject to inference. In other words, if the user does not specify any mock model parameters, the generative mock model will be the same as the model used to fit the data. In most cases, one will be interested in studying mismodeling, i.e., when the mock data is generated from a model different from that used to fit the data. This can be achieved by deliberately setting the mock model parameters different from the prior model. 

Selecting the initial state
+++++++++++++++++++++++++++++
When the dataset is supplied by the user (``strgexprflux`` is set), and unless specified otherwise by setting ``randinit=False``, the initial state of the chain is drawn randomly from the prior. Note that in order for the initial state to be nonrandom, a reference catalog needs to be internally supplied.

In constrast, if the dataset is simulated, ``randinit`` is not ``True`` and generative mock model is the same as the prior model, then the initial state of the chain is set to be the state from which the mock dataset was drawn.  


.. _sectoutp:

Output
-------------
A function call to ``pcat.main.init()`` returns the collected samples as well as postprocessed variables in an object that we will refer to as ``gdat``. Any output (as well as many internal variables of the sampler) can be accessed via the attributes of this global object. 

Furthermore, PCAT ships with extensive routines to visualize the output chain. The output plots are placed in the relevant subfolders ``pathbase/data/outp/rtag`` and ``pathbase/imag/rtag``, respectively, where ``rtag`` is the run tag. The ``pathbase`` folder is created if it does not already exist. It defaults to the value of the environment variable ``$PCAT_DATA_PATH``. Therefore ``pathbase`` can be omitted by setting the environment variable ``$PCAT_DATA_PATH``.


.. _sectplot:

Plots
+++++
If not disabled by the user, PCAT produces plots in every stage of a run. Some plots are produced in the initial setup, frame plots are produced at predetermined times during the sampling and others are produced in the postprocessing after all chains have run. There are five subfolders in the plot folder of a given run.

- ``init`` Initial setup, pixel lookup table (if applicable)
- ``fram`` Frame plots
- ``diag`` Diagnostic plots
- ``post`` Posterior distribution plots of model parameters and derived quantities, prior and likelihood.
- ``anim`` GIF animations made from the frame plots in ``fram`` that are produced during sampling.


Chain
+++++
``pcat.main.init()`` returns an object that contains the output chain.


Diagnostics
------------

Autocorrelation
+++++++++++++++++++++++++++++

A chain of states needs to be Markovian (memoryless) in order to be interpreted as fair draws from a target probability density. The autocorrelation of the chain shows whether the chain is self-similar along the simulation time (either due to low acceptance rate or small step size). Therefore the autocorrelation plots should be monitored after eah run.

In a transdimensional setting, the autocorrelation of a parameter is ill-defined, since parameters can be born, killed or change identity. Therefore, for such parameters, we calculate the autocorrelation
The autocorrelation Note that the
In the The acceptance rate of the birth and death moves shows whether the prior on the amplitude of the elements is appropriate. If the acceptance rate of birth and deaths proposals is too low, the lower limit of the prior on the amplitude is too high. Likewise, if it is too low, the lower limit of the prior on the amplitude is too low. This behaviour is due to the fact that element parameters (position, amplitude and, if relevant, spectral or shape parameters) are drawn randomly from the prior when a birth is proposed. Therefore the lower limit on the amplitude prior should be adjusted such that the birth and death acceptance rate is between :math:`\sim 5%% - \sim30%%`. Otherwise across-model sampling will be inefficient and result in slow convergence.  

Gelman-Rubin test
+++++++++++++++++++++++++++++
PCAT nominally runs multiple, noninteracting chains, whose samples are aggregated at the end. In order to ensure convergence, therefore, one can compare within-chain variance with across-chain variance. This is known as the Gelman-Rubin test. PCAT outputs the GR test statistics in ``gdat.gmrb`` and plots the relevant diagnostics in ``$PCAT_DATA_PATH/imag/rtag/diag/``.

.. note:: Make sure to run PCAT with the argument ``diagmode=False`` for reasonable time performance. ``diagmode=True`` option puts the sampler in a conservative diagnostic mode and performs extensive checks on the state of critical data structures to ensure that the model and proposals are self consistent, largely slowing down execution.


Tutorial
--------------
In this tutorial, we will illustrate how to run PCAT during a typical science analysis. Assuming that you have :ref:`installed <sectinst>` PCAT, let us first run it on mock data.

All user interaction with PCAT can be performed through the ``pcat.main.init()`` funtion. Because arguments to this function have hierarchically defined defaults, even the default call (without arguments) starts a valid PCAT run.

.. code-block:: python

    import pcat.main
    pcat.main.init()

This runs PCAT on mock Fermi-LAT data to collect a single chain of 100000 samples before thinning and burn-in. After initialization, PCAT collects samples, produces frame plots (snapshots of the sampler state during the execution) and postprocesses the samples at the end. The run should finish in under half an hour with the message

.. code-block:: none
    
    >> The ensemble of catalogs is at $PCAT_DATA_PATH/data/outp/rtag/

While the sampler is running, you can check ``$PCAT_DATA_PATH/imag/rtag/`` to inspect :ref:`the output plots <sectplot>`.


.. code-block:: python

    pcat.main.init( \
              randinit=True, \
              maxmgang=deg2rad(20.), \
              indxenerincl=arange(1, 4), \
              indxevttincl=arange(2, 4), \
              lgalcntr=lgalcntr, \
              bgalcntr=bgalcntr, \
              minmflux=3e-11, \
              maxmflux=3e-6, \
              strgback=['isotflux.fits', 'fdfmflux%s.fits' % strgcntr], \
              strgexpo='fermexpo_cmp0_igal%s.fits' % strgcntr, \
              strgexprflux='fermflux_cmp0_igal%s.fits' % strgcntr, \
             )
    
    

API
---------

All user interaction with PCAT is accomplished through the ``pcat.main.init()`` function. Below is a list of its function arguments.

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


    **Associations with the reference catalog**
 
    :param anglassc: Radius of the circle within which sample catalog point sources can be associated with the point sources in the reference catalog.
    
    :type anglassc: float


    :param margfactcomp: The ratio of the side of the square in which the sample point sources are associated with the reference catalog, to the size of the image.

    :type margfactcomp: float


    :param nameexpr: A string that describes the provided reference catalog to be shown in the plot legends.

    :type nameexpr: str
    

    :param cntrpnts: Force the mock data to a single PS at the center of the image. Defaults to ``False``.
    
    :type cntrpnts: bool


    **General**

    :param pntstype: Functional type of point sources. 
        
        - ``lght`` Point sources are light sources
        - ``lens`` Point sources are lenses.
    
    :type pntstype: str


    **Proposals**

    :param propcova: Propose from the full covariance, taking steps in all parameters in a single proposal. Otherwise, a step is taken along only one dimension. Defults to False.

    :type randinit: bool


    **Initial state**

    :param randinit: Force the initial state to be randomly drawn from the prior. Default behavior for mock data is to initialize the chain with the true state.

    :type randinit: bool


    **Adaptive burn-in**

    :param loadvaripara: Load the diagonal elements of the previously learned covariance

    :type loadvaripara: bool


    :param optiprop: Optimize the scale of each proposal by acceptance rate feedback. All samples during the tuning are discarded.
        
    :type optiprop: bool


    **Post processing**

    :param regulevi: Regularize the log-evidence estimate.

    :type regulevi: bool


    :param strgexprflux: Name of the FITS file (without the extension) in ``pathdata`` containing the observed data as an ``ndarray``.

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


    :param asscmetrtype: Type of metric used to associate the sample catalogs with the reference catalog

    :type asscmetrtype: str


    **Diagnostics**

    :param emptsamp: Perform a futile run without collecting any samples, but creating all data structures and producing all visualizations as if in a normal run. Defaults to ``False``.

    :type emptsamp: bool


    :param diagmode: Start the run in diagnostic mode. Defaults to ``False``.

    :type diagmode: bool


    **Plotting**

    :param strgexprname: A string describing the experiment used to collect the observed data.

    :type strgexprname: str


    :param strganglunit: Label for the spatial axes.

    :type strganglunit: str

    
    :param labllgal: Label for the horizontal axis

    :type labllgal: str

    
    :param lablbgal: Label for the vertical axis

    :type lablbgal: str

    
    :param spatdisttype: Type of spatial distribution of sources for each population

    :type spatdisttype: list of str
    

    :param fluxdisttype: Type of flux distribution of sources for each population

    :type fluxdisttype: list of str
    

    :param spatdisttype: Type of spatial distribution of sources for each population

    :type spatdisttype: list of str
    

    :param spectype: Type of energy spectrum of sources for each population

    :type spectype: list of str
    

    :param psfntype: Type of PSF radial profile

    :type spectype: list of str
    

    :param oaxitype: Type of PSF off-axis profile

    :type oaxitype: list of str
    

.. ## PSF
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
         psfninfoprio=True, \
         ## spectral
         # prior
         priotype='logt', \
         priofactdoff=0., \
         margfactmodl=0.9, \
         bindprio=False, \
         
         maxmgang=None, \
         minmmeanpnts=None, \
         gangdistscal=None, \
         bgaldistscal=None, \
         curvdistmean=None, \
         curvdiststdv=None, \
         minmflux=None, \
         minmfluxdistslop=None, \
         minmfluxbrek=None, \
         minmfluxdistbrek=None, \
         minmfluxdistsloplowr=None, \
         minmfluxdistslopuppr=None, \
         minmsinddistmean=None, \
         minmsinddiststdv=None, \
         
         minmsigm=None, \
         meansigm=None, \
         stdvsigm=None, \
         minmgamm=None, \
         meangamm=None, \
         stdvgamm=None, \
         minmpsff=None, \
         meanpsff=None, \
         stdvpsff=None, \
         
         bacp=None, \
         fluxsour=None, \
         fluxsour=None, \
         sizesour=None, \
         sizesour=None, \
         ratisour=None, \
         ratisour=None, \
         ellphost=None, \
         ellphost=None, \
         sherhost=None, \
         sherhost=None, \
         beinhost=None, \
         beinhost=None, \
         
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
         numbsidecart=200, \
         numbsideheal=256, \
         numbdatasamp=100, \

    **Unfunctional**
 
    :param numbspatdims: Number of spatial dimensions. Currently not functional.

    :type numbspatdims: None


.. note::
    
    The mock model parameters can be set by preceeding the parameter name with ``mock``. For example, in order to set the mock number of PS, you can specify ``mocknumbpnts=array([10])``.



    


Garbage collection
-------------------

PCAT produces two folders for the output of each run, one for plots and the other to contain the chain saved to the disc. Given that many test and intermediate runs may be needed before each science run, the number of folders (and files therein) may increase quickly. In order to avoid this, the script ``gcol.py`` is provided to the user as a convenience. When executed, it erases the output from all runs that

- have not run to completion, i.e., does not have the animations, which are produced at the very end, or
- has collected less than 100000 samples per chain.

