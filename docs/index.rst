Welcome to PCAT's documentation!
================================

When testing hypotheses or inferring their free parameters, a recurring problem is to compare models that contain a number of objects whose multiplicity is itself unknown. Therefore, given some data, it is desirable to be able to compare models with different number of parameters. One way of achieving this is to obtain a point estimate (usually the most likely point) in the parameter space of each model and, then, rely on some information criterion to penalize more complex models for excess degrees of freedom. Another way is to sample from the parameter space of each model and compare their Bayesian evidence. Yet another is to take samples from the union of the hypotheses using a set of transdimensional jumps across models. This is what PCAT (Probabilistic Cataloger) is designed for.

PCAT is a hierarchical, transdimensional MCMC sampler. Given some data, it can be used to sample from the **catalog space**, i.e., the hypothesis space of a model with a variable number of components. Alternatively, it can be used as a general purpose Poisson mixture sampler to infer clusters in a dataset.

A user manual for PCAT will be published here soon.



.. toctree::
    :maxdepth: 2
    :caption: Installation:
    Welcome to PCAT's documentation
    installation

Installation
==================

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
