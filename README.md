# PCAT (Probabilistic Cataloger)

PCAT is a Bayesian framework to sample from the catalog space. 

## Usage

### Options
---
#### Diagnostics
`diagsamp`
Boolean flag to run the sampler in diagnostic mode

---
#### User interaction
`verbtype`
Verbosity level
- `0`: no standard output
- `1`: print only critical statements along with periodic progress
- `2`: full log of internal variables, to be used only for diagnostic purposes

---
#### Plotting
`numbswepplot`
Number of sweeps between frames

`makeplot`
Boolean flag to allow making of plots

`trueinfo`

---
#### Sampler
`numbproc`

`numbswep`
Number of sweeps

`numbburn`
Number of sweeps to be discarded during burn-in
         
`factthin`
The factor by which to thin the chain

`probprop`

`datatype`
Type of data
- `'inpt'`: provided by the user
- `'mock'`: generated (mock) data

`regitype`
Type of region
- `'igal'`: ROI centered at the galactic center
- `'ngal'`: ROI centered at the North Galactic Pole


---
#### Initial state
`initnumbpnts`

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

`maxmnormback`

`minmnormback`

---
#### Hyperpriors
`minmfdfnnorm`

`maxmfdfnnorm`

`fdfntype`

`minmfdfnslop`

`maxmfdfnslop`

`minmfdfnbrek`

`maxmfdfnbrek`

`minmfdfnsloplowr`

`maxmfdfnsloplowr`

`minmfdfnslopuppr`

`maxmfdfnslopuppr`

`meansdfn`

`stdvsdfn`

---
#### PSF
`psfntype`

`proppsfn`

`maxmangleval`

---
#### General
`maxmnumbpnts`

`liketype`

`exprtype`

---
#### ROI
`pixltype`

`lgalcntr`

`bgalcntr`

`margsize`

`numbsidecart`

`numbsideheal`

---
#### Proposal scales
`stdvfdfnnorm`

`stdvfdfnslop`

`stdvpsfipara`

`stdvback`

`stdvlbhl`

`stdvflux`

`stdvsind`

`fracrand`

`radispmrlbhl`

---
#### Mock data
`mocknumbpnts`

`mockfdfnslop`

`mockfdfnsloplowr`

`mockfdfnslopuppr`

`mockfdfnbrek`

`mocknormback`

`mockfdfntype`

`mockpsfntype`

---
#### Exposure
`strgexpo`

---
#### User-provided aata
`strgexpr`
`indxevttincl`

`indxenerincl`

---
#### Background modeling
`strgback`

`lablback`

`nameback`

