# PCAT (Probabilistic Cataloger)

PCAT is a Bayesian framework to sample from the catalog space. It's theoretical framework is introduced in [Daylan, Portillo & Finkbeiner (2016)](https://arxiv.org/abs/1607.04637), submitted to ApJ. Refer to its [webpage](http://www.tansudaylan.com/pcat) for an introduction.


### Installation

You can install `tdpy` by running the `setup.py` script.
```
python setup.py install
```

### Usage
The basic interaction with PCAT is through the `init` function, which writes the resulting probabilistic catalog to the disc. The software ships with convenience functions to post-process the output.

```python
import pcat
     
numbener = 5
minmspec = array([3e-11])
maxmspec = array([1e-7])
mockfluxdistslop = array([[1.9]])
pcat.init(, \
          psfntype='doubking', \
          randinit=False, \
          trueinfo=True, \
          maxmgang=20., \
          mocknumbpnts=array([300]), \
          minmspec=minmspec, \
          maxmspec=maxmspec, \
          maxmnormback=array([2., 2.]), \
          minmnormback=array([0.5, 0.5]), \
          strgexpo='fermexpo_comp_ngal.fits', \
          datatype='mock', \
          numbsideheal=256, \
          mockfluxdistslop=mockfluxdistslop, \
          mocknormback=ones((2, numbener)), \
         )
```


### Options
---
#### Diagnostics
`diagmode`
Boolean flag to run the sampler in diagnostic mode

---
#### User interaction
`verbtype`
Verbosity level

---
#### Plotting
`numbswepplot`
Number of sweeps between frames

`makeplot`
Boolean flag to allow making of plots

`trueinfo`
Boolean flag to allow cross-correlation against the true data and superposition of true values in the sample and posterior plots

---
#### Sampler
`numbproc`
Number of processes

`numbswep`
Number of sweeps

`numbburn`
Number of sweeps to be discarded during burn-in
         
`factthin`
The factor by which to thin the chain

`probprop`
Array indicating the probability of proposing the associated proposal types

`datatype`
Type of data

---
#### Initial state
`initnumbpnts`
Initial number of Point Sources (PSs)

`randinit`
Boolean flag to start the MCMC state randomly from the prior

`maxmgang`
Half-size of the ROI window

`minmflux`
Minimum PS flux allowed by the model

`maxmflux`
Maximum PS flux allowed by the model

`minmsind`
Minimum spectral index allowed by the model

`maxmsind`
Maximum spectral index allowed by the model

`minmnormback`
Minimum background normalizations

`maxmnormback`
Maximum background normalizations

---
#### Hyperpriors
`minmmeanpnts`
Minimum mean number of PS

`maxmmeanpnts`
Maximum mean number of PS

`fluxdisttype`
Flux Distribution Function (FDF) shape

`minmfluxdistslop`
Minimum FDF slope

`maxmfluxdistslop`
Maximum FDF slope

`minmfluxdistbrek`
Minimum FDF break flux

`maxmfluxdistbrek`
Maximum FDF break flux

`minmfluxdistsloplowr`
Minimum FDF lower slope

`maxmfluxdistsloplowr`
Maximum FDF lower slope

`minmfluxdistslopuppr`
Minimum FDF upper slope

`maxmfluxdistslopuppr`
Maximum FDF upper slope

`sinddistmean`
Mean Spectral index Distribution Function (SDF)

`sinddiststdv`
Standard deviation of the SDF

---
#### PSF
`psfntype`
Functional form of the PSF

`boolproppsfn`
Boolean flag to propose PSF updates

`boolpropfluxdist`
Boolean flag to propose FDF updates

`maxmangleval`
Maximum angular distance up to which point source flux will be calculated.

---
#### General
`maxmnumbpnts`
Maximum number of point sources allowed

`liketype`
Likelihood to be evaluated

`exprtype`
Experiment type, which allows experiment-specific defaults to be loaded

---
#### ROI
`pixltype`
Pixelizaiton type

`lgalcntr`
Longitude of the ROI center

`bgalcntr`
Latitude of the ROI center

`margsize`
Size of the margin around the ROI, where model point sources can exist

---
#### Proposal scales
`stdvmeanpnts`
Proposal scale for the mean number of point sources

`stdvfluxdistslop`
Proposal scale for the slope of the FDF

`stdvpsfp`
Proposal scale for the PSF parameters

`stdvback`
Proposal scale for the mean number of point sources

`stdvlbhlminm`
Minimum proposal scale for spatial coordinates

`stdvlbhlmaxm`
Maximum proposal scale for spatial coordinates

`stdvflux`
Proposal scale for fluxes

`stdvsind`
Proposal scale for spectral indices

`fracrand`
Probability of proposing the next state randomly from the prior

`radispmr`
Radius of the circle in which splits and merges will be proposed

---
#### Mock data
`mocknumbpnts`
Mock number of point sources

`mockfluxdistslop`
Mock FDF slope

`mockfluxdistsloplowr`
Mock FDFlower slope

`mockfluxdistslopuppr`
Mock FDF upper slope

`mockfluxdistbrek`
Mock FDF break flux

`mocknormback`
Mock background normalization

`mockfluxdisttype`
Mock FDF shape

`numbsidecart`
Number of pixels along the side of the ROI for the generated data map, if the pixelization is Cartesian.

`numbsideheal`
Base resolution of the HealPix map, if the pixelization is HealPix.

---
#### Exposure
`strgexpo`
File name of the backgrounds in the data path

---
#### User-provided aata
`strgexpr`
File name of the expsoure in the data path

`indxevttincl`
Indices of the PSF classes to be included in the fit

`indxenerincl`
Indices of the energy bins to be included in the fit

---
#### Background modeling
`strgback`
File name of the backgrounds in the data path

`lablback`
Mathematical label of the backgrounds

`nameback`
Descriptive name of the backgrounds
