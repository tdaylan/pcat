# PCAT (Probabilistic Cataloger)

## Installation

PCAT is a transdimensional, hierarchical, and Bayesian framework to sample from the posterior probability distribution of a metamodel (union of models with different dimensionality) given some Poisson-distributed data. 

In addition to its previous use in the literature to sample from [the point source catalog space consistent with the Fermi-LAT gamma-ray data](http://iopscience.iop.org/article/10.3847/1538-4357/aa679e/meta), and [the dark matter subhalo catalog space consistent with an Hubble Space Telescope (HST) optical image](https://arxiv.org/abs/1706.06111), it can also be used as a general-purpose Poisson mixture sampler.

During burn-in, it adaptively optimizes its within-model proposal scale to minimize the autocorrelation time. Furthermore, it achieves parallelism through bypassing Python's Global Interpreter Lock (GIL). It is implemented in ```python2.7``` and its theoretical framework is introduced in [Daylan, Portillo & Finkbeiner (2016)](https://arxiv.org/abs/1607.04637). Refer to its [webpage](http://www.tansudaylan.com/pcat) for an introduction.


## Installation

You can install PCAT either by using ```pip```
```
pip install pcat
```

or, by running its `setup.py` script.

```
python setup.py install
```

Note that PCAT depends on [TDPY](https://github.com/tdaylan/tdpy), a library of MCMC and numerical routines. The pip installation will install PCAT along with its dependencies.

## Usage
PCAT user manual is on [ReadTheDocs](http://pcat.readthedocs.io/en/latest/).

