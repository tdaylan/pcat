Likelihood configurations
=========================

PCAT separates data likelihoods from transdimensional catalog postprocessing.
Every transdimensional model produces a catalog at each retained draw. With
the default ``boolcondcatl=True``, final processing builds a condensed posterior
catalog by associating recurrent supported spatial elements across those
sampled catalogs. Optional matching to an external reference catalog follows
this summary. Condensation and either kind of association are posterior
products, not likelihoods. Set ``boolcondcatl=False`` to skip condensation.
The current implementation does not spatially condense one-dimensional line,
flare, or Keplerian catalogs.

Poisson image likelihood
------------------------

The image workflow predicts counts in pixels from a source model, exposure,
background emission, and instrument response. Point sources, extended emission,
foreground light, and gravitational lenses are components of this one image
likelihood, not separate likelihood families.

Data and pixelization
~~~~~~~~~~~~~~~~~~~~~

Set ``typedata="simu"`` for generated counts or ``typedata="inpt"`` for
supplied data. For input data, ``strgexprsbrt`` identifies the surface-brightness
or count-rate FITS input. Exposure can be constant with ``typeexpo="cons"``
and ``strgexpo`` or read from a FITS map with ``typeexpo="file"`` and
``strgexpo``. PCAT converts the input to expected counts using the exposure and
its configured units before evaluating the likelihood.

``typepixl="cart"`` selects a Cartesian image. Set ``numbsidecart`` to choose
its side length; ``boolforccart=True`` forces a square grid when pixel
coordinates are supplied directly. ``typepixl="heal"`` selects HEALPix data
and uses ``numbsideheal`` for its resolution. Fermi-LAT defaults to HEALPix;
other image experiments default to Cartesian pixels.

``liketype="pois"`` is the default for count data. Omitting the data-only
factorial, its log likelihood is

.. math::

   \log L_{\mathrm{Pois}} = \sum_i D_i \log M_i - M_i,

where ``D`` is the observed count map and ``M`` is the predicted count map.
Set ``liketype="gaus"`` for the built-in Gaussian residual likelihood, which
uses ``varidata`` (by default, the observed counts with a floor of one) as its
per-pixel variance:

.. math::

   \log L_{\mathrm{Gaus}} = -\frac{1}{2}\sum_i \frac{(D_i-M_i)^2}{V_i}.

Source and lens components
~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure populations with ``typeelem`` and their true and fitted settings in
``dicttrue`` and ``dictfitt``. Common image components include:

* ``lghtpnts`` for point sources with sampled positions and fluxes.
* ``lghtgausbgrd`` for extended Gaussian emission components.
* ``lens`` for lens deflectors and their associated mass parameters.
* ``clusvari`` for variable-width Gaussian components represented in an image count map.
* ``typeemishost="sers"`` for a Sérsic foreground lens-galaxy light model;
  ``"none"`` disables this host-emission component.
* ``boollenssubh=True`` to include lens-substructure deflectors where the
  selected lens configuration supports them.

The ``typespatdist`` or ``spatdisttype`` setting selects the spatial prior for
each population. Source flux priors are selected with ``typeprioflux`` and
bounded or customized with the corresponding per-element limits and
``limtparaelem`` where available. HST strong-lens examples combine a fixed
macro lens and source/host light with a variable ``lens`` population for
additional deflectors.

Background emission and exposure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Experiment defaults define the available background components. Override their
templates with ``sbrtbacknorm`` in ``dicttrue`` or ``dictfitt``. The
``listnamediff`` entry selects the diffuse components added to the model; names
such as ``back0000`` refer to entries in the template list. ``boolspecback``
controls whether a background amplitude is shared spectrally or varies by
energy bin. Fermi-LAT configurations can use a uniform component, a diffuse
template, or the ``anlytype="bfun"`` background-function basis. Exposure
controls how model brightness becomes expected counts and can be masked with
``listmask`` and ``typemaskexpo``.

Point-spread functions
~~~~~~~~~~~~~~~~~~~~~~

``typeevalpsfn`` controls where the PSF is applied:

* ``"none"`` leaves the source model unconvolved.
* ``"kern"`` applies a local PSF kernel to supported point-source elements.
* ``"conv"`` convolves eligible image components.
* ``"full"`` enables both kernel evaluation and image convolution where
  supported.

The experiment selects a PSF profile family through ``typemodlpsfn``. Current
experiment defaults include single Gaussian, single King, and double King
profiles. Supply profile parameters through ``psfpexpr`` in ``dicttrue`` and
``dictfitt``; set ``boolmodipsfn=True`` when PSF parameters should vary in the
fit. Fermi-LAT can also use its energy- and event-class-dependent PSF scale
``fermscalfact``. PSF and source-kernel options depend on pixelization and
element type; unsupported combinations should use ``typeevalpsfn="none"``.

Unbinned event-level likelihoods
--------------------------------

PCAT can sample an unbinned likelihood through a user callback. The maintained
``unbinned_gaussian_mixture`` example passes the 140 simulated event coordinates
directly to ``pcat.fixed.sample_fixed_chains`` and sums their Gaussian-mixture
log densities with no spatial bins or count map. It samples the four center
coordinates of two equal-weight Gaussian components; the component count and
width are fixed. This is a fixed-dimensional custom-likelihood demonstration,
not a native transdimensional Gaussian-mixture likelihood.

The native ``typeexpr="gmix"`` experiment remains the binned Poisson workflow
shown in ``gaussian_mixture_catalog``. Its variable-width ``clusvari`` elements
are fit to an 8-by-8 count grid. ``boolbinsener=False`` only indicates that this
example has no energy axis. PCAT does not currently provide a built-in unbinned
``typeexpr="gmix"`` dispatcher for variable-count mixture components.

Spectral and time-series likelihood
-----------------------------------

The one-dimensional ``typeexpr="fire"`` workflow uses ``anlytype="spec"`` for
spectral bins. Set ``binsenerfull`` to the axis edges and use
``typedata="inpt"``, ``strgexprsbrt``, ``typeexpo``, and ``strgexpo`` for
supplied spectra and exposure. A continuum or background template can be
provided through ``dictfitt["sbrtbacknorm"]`` and activated with
``dictfitt["listnamediff"]``.

The ``spectype`` setting selects an analytic spectral shape:

* ``"powr"`` for a power law.
* ``"gaus"`` for a Gaussian profile.
* ``"lore"`` for a Lorentzian profile.
* ``"voig"`` for a Voigt profile.
* ``"pvoi"`` for a pseudo-Voigt profile with mixing fraction ``frac``.
* ``"sinc"`` for a sinc-squared profile.
* ``"skew"`` for a skew-Gaussian profile with shape parameter ``skew``.
* ``"toph"`` for a top-hat profile.
* ``"edis"`` for an energy-dispersion profile using the configured
  ``edisintp`` response.
* ``"colr"``, ``"curv"``, or ``"expc"`` for color-index, curved, or
  exponential-cutoff continua.

Line elements use integrated flux ``flux`` and line center ``elin``. Gaussian,
sinc-squared, skew-Gaussian, and top-hat profiles add ``sigm``. Lorentzian
profiles add ``gamm``. Voigt profiles use both widths. Pseudo-Voigt profiles
also use ``frac``. Set parameter bounds with ``limtparaelem`` and proposal
scales with ``stdvpropelemfire``.

Set ``lsftype="gaus"`` and ``lsfresolvingpower`` to convolve each line at a
constant instrument resolving power. Set ``lsftype="tabu"`` and pass an
odd-length, nonnegative ``lsfkernel`` for a measured line-spread function.
Both paths preserve sampled line flux. The default ``lsftype="none"`` leaves
the intrinsic profiles unchanged.

The same one-dimensional machinery can model time-binned data by treating time
as the axis. ``flargauss`` uses peak amplitude ``flux``, peak time ``elin``,
and full width at half maximum ``fwhm``. ``flarexpd`` adds a one-sided decay
time ``scalfall``. ``flarfred`` uses independent ``scalrise`` and ``scalfall``
for a fast-rise, exponential-decay profile. ``flardav`` uses ``fwhm`` for the
empirical Davenport flare template. The maintained stellar-flare example fits
native ``flarfred`` elements. Radial-velocity models use Keplerian
``lghtlinekepl`` elements and their dedicated likelihood in
:mod:`pcat.radial_velocity`.

Transit-timing likelihood
-------------------------

For transit-time data, use ``pcat.time_series.log_likelihood_transit_times``
inside a fixed-dimensional PCAT likelihood callback. Pass observed transit
times, predicted transit times from Ephesos, and positive timing errors in the
same units. The function returns the uncertainty-weighted Gaussian log
likelihood up to a parameter-independent normalization constant; it returns
negative infinity for nonfinite model times. Ephesos supplies linear,
sinusoidal-residual, and optional N-body time predictions. PCAT handles the
parameter priors and posterior sampling through ``sample_fixed`` or
``sample_fixed_chains``. The predictor must align observed transit epochs with
model epochs before evaluating this likelihood.

User-supplied likelihood
------------------------

Use ``typeexpr="gener"`` for a fixed-dimensional custom likelihood. The
``retr_llik(gdat, strgmodl, values)`` callback returns one log-likelihood value
for the supplied parameter vector; PCAT handles configured priors, proposals,
accept/reject steps, persistence, and diagnostics. This is the appropriate path
when none of the built-in binned likelihoods matches the data model. See
:doc:`capabilities` for a complete callback example and :doc:`sampling` for
proposal and convergence controls.