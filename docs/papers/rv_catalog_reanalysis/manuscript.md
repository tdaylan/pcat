# Transdimensional Reanalysis of Public Exoplanet Radial-Velocity Curves

**Working single-author manuscript draft.** The data inventory and initial pipeline pilots are complete. The full per-curve posterior analysis and literature-reference matching are not. Do not interpret this draft as a completed scientific result.

Tansu Daylan  
[Affiliation to confirm]

## Abstract

Published stellar radial-velocity (RV) curves are commonly interpreted through one or more fixed-dimensional Keplerian solutions, although the number of signals supported by a given time series can be uncertain. We develop a reanalysis of public RV curves using PCAT, a Bayesian sampler that treats the number and parameters of Keplerian signals as a transdimensional catalog. The reproducible input boundary is the three RV time-series collections listed by the NASA Exoplanet Archive on 2026-10-08: the literature RADIAL manifest, MARVELS, and the Howard--Fulton completeness sample. The archive lists 1,072 RADIAL manifest entries, with one repeated filename; the snapshot contains 1,071 unique curves for 511 archive identifiers, 6,320 MARVELS curve files for 3,169 header target identifiers, and 76 Fulton stellar files. A file-level audit finds inconsistent time-reference labels, repeated target/epoch keys, conflicting reductions, and non-finite MARVELS rows. We therefore preserve source curves as separate likelihood units rather than silently combining them. Run-length pilots on the 53-observation HD 142 curve give materially different planet-count posteriors, and none establishes convergence or supports revising the published system solution. The paper specifies a reproducible archive curation and per-curve analysis protocol; full posterior results, sensitivity tests, and comparison against the source publications remain necessary before scientific conclusions can be drawn.

## 1. Introduction

Stellar radial velocities constrain the gravitational reflex motion of a host star and remain a primary route to discovering and characterizing exoplanets. The inferential object is not always a fixed list of planets. Signals can be weak, periods can be aliased by the observing window, eccentricity can covary with phase and amplitude, and offsets between instruments can compete with long-period trends. A conventional fit selects a planet count and estimates orbital parameters conditional on that count. If several counts explain the observations, a single best-fit solution does not represent that uncertainty.

PCAT applies probabilistic cataloging to this setting by representing each Keplerian signal as an exchangeable catalog element. The sampler explores models with different numbers of planets and orbital parameters in a single transdimensional chain. The RV likelihood analytically marginalizes linear instrument offsets and can numerically marginalize stellar jitter. Birth and prior-redraw jump proposals use a periodogram-weighted density with the matching Hastings correction, improving access to period regions suggested by the data without making a periodogram peak a detection claim.

This study asks a limited, testable question: for each public archive RV curve, how does the posterior over the number and orbital parameters of Keplerian signals compare with the interpretation in the source publication? We do not treat the archive as a volume-limited or detection-complete planet sample. Its constituent data sets are literature- and survey-selected, and some targets occur in multiple reductions. The unit of initial inference is therefore the archived curve, not an assumed independent planet system.

## 2. Public RV data and curation

### 2.1 Archive snapshot

The NASA Exoplanet Archive's bulk-download listing dated 2026-10-08 identifies three relevant public products. RADIAL is a set of literature RV curves, MARVELS is a separate survey release, and Fulton provides RV observations used in a survey-completeness study. As retrieved on 2026-10-11, these products contain 1,071 unique RADIAL curve files for 511 archive UIDs, 6,320 MARVELS files with 3,169 distinct header `STAR_ID` values, and 76 Fulton stellar series. The RADIAL manifest itself has 1,072 download lines because one filename is listed twice.

| Collection | Curve files | Target identifiers | Parsed RV rows | Valid rows after finite-value and positive-error checks |
|---|---:|---:|---:|---:|
| RADIAL | 1,071 unique (1,072 manifest entries) | 511 archive UIDs; 552 header name strings | 43,704 | 43,474 |
| MARVELS | 6,320 | 3,169 header `STAR_ID` values | 171,351 parsed; one header-point discrepancy | 138,076 |
| Howard--Fulton | 76 | 76 HD identifiers | 9,613 | 9,613 |

The RADIAL curves use ten reported time-reference labels and 68 distinct instrument strings. Among valid RADIAL records, 1,444 target/time-frame/epoch keys recur across files. At the recorded precision, 225 row pairs have identical velocities and uncertainties, while 1,344 pairs disagree in at least one of those quantities. MARVELS contains 33,275 rows excluded by the finite-value and positive-uncertainty rule, plus one unparsable row or header count discrepancy. These counts are an audit of this downloaded snapshot, not counts of independent observations or stars in a uniform survey.

### 2.2 Inclusion and exclusions

The preparation script stores the raw collection, original curve filename, target identifiers, source reference or Bibcode, telescope, instrument, time-reference label, values, uncertainties, and a row-validity reason. It does not rewrite source values, infer missing uncertainties, resolve aliases, convert time standards, or merge cross-file observations. Fulton files retain their row-level instrument codes. All raw downloads and normalized tables remain outside the source repository; the script records the RADIAL manifest hash and source-file hashes.

Under the current runner's minimum of ten valid observations and same-curve repeated-epoch rule, 7,051 of the 7,467 files are eligible for fitting, 413 have too few valid rows, and three contain conflicting measurements at a repeated epoch. Eligibility does not imply that a curve's time system, target identity, or noise model is scientifically resolved.

Each curve is analyzed independently in the initial pass. This prevents conflicting published reductions from being multiplied as if they were independent measurements. Repeated target/epoch keys across curves are reported for later curation. No aggregate multiplicity frequency is calculated from these data because archive inclusion and observation selection are not a defined survey selection function.

## 3. PCAT radial-velocity model

For a time series from one archive curve, the signal model is a sum of $N$ Keplerian reflex-velocity components,

$$
v(t_i) = \sum_{n=1}^{N} v_n(t_i; P_n, K_n, e_n, \omega_n, M_{0,n}) + \sum_{j=1}^{J} c_j X_{ij} + \epsilon_i.
$$

Here $P_n$ is period, $K_n$ is semi-amplitude, $e_n$ is eccentricity, $\omega_n$ is argument of periastron, and $M_{0,n}$ is mean anomaly at the reference epoch. The design matrix $X$ contains one constant offset per instrument and, when configured, a linear trend. These linear nuisance coefficients are marginalized analytically under flat priors. An additional white-noise jitter term is marginalized numerically over the RV adapter's log-spaced jitter grid. The row uncertainty supplied by the archive remains the reported measurement uncertainty.

PCAT samples $N$ between the configured lower and upper bounds using birth and death transitions, and updates existing orbital elements with within-model and prior-redraw jump proposals. The RV adapter also provides a data-informed birth/jump proposal built from a weighted periodogram-like search over circular fits. The proposal density, including its uniform mixture component, is included in the Hastings ratio. Data-informed proposals guide exploration; they are not used as a detection statistic.

The initial full-analysis configuration permits zero through five Keplerian elements, uses the archive curve's native time values within each curve, supplies one instrumental label per reported instrument, and enables likelihood-tempered burn-in followed by ordinary target sampling. Prior bounds, jitter grid, proposal scales, random seeds, software revision, and run length will be recorded for every curve. Curves with insufficient finite measurements or unresolved repeated epochs within the same file are retained in the audit table and assigned an explicit non-fit status rather than silently altered.

## 4. Analysis and interpretation plan

For each fitted curve, the primary result is the posterior probability of signal count $P(N\mid d)$ for $0\leq N\leq5$. Secondary summaries are posterior distributions of period and semi-amplitude, posterior-predictive RV curves, likelihood and model-count traces, proposal acceptance, and convergence/effective-sample diagnostics. Source references in the archive metadata will be matched to the corresponding reported planet solution when the target and reference can be resolved. Agreement is assessed at the curve/reference level; repeated curves for one target are not pooled until overlap and calibration are resolved.

A result is considered informative only if the retained chain explores the supported model-count range relevant to that curve, posterior summaries are stable to longer runs or independent seeds, and the result is not dominated by prior bounds. Curves failing these conditions remain unresolved. We will report the number of curves downloaded, parseable, eligible, completed, converged under predeclared criteria, and unresolved, by archive collection. The per-curve output table and exact inputs will accompany the manuscript.

### 4.1 Pilot check

One RADIAL curve for HD 142 (archive file `UID_0000522_RVC_001.tbl`, 53 valid observations, source Bibcode `2006ApJ...646..505B`) was used for initial run-length pilots. A 2,000-sweep run with 1,000 retained draws assigns 0.93 to $N=2$ and 0.07 to $N=3$. Other 2,000-sweep pilot configurations yielded $P(N=1,2,3)=(0.39, 0.54, 0.07)$ or $P(N=2,3)=(0.93, 0.07)$. Two 100,000-sweep runs with independent seeds and 2,000 retained draws give $P(N=2,3,4,5)=(0.753, 0.2165, 0.0295, 0.001)$ and $(0.840, 0.1415, 0.0175, 0.001)$. Their count occupancy changes across chain halves; for seed 17, $P(N=2)$ is 0.863 in the first half and 0.643 in the second, while seed 18 gives 0.811 and 0.869. Thus the long-run chains broadly agree on a two-planet-dominated posterior but do not establish convergence. These pilots do not support revising the published HD 142 solution.

## 5. Limitations

This archive snapshot is not “all public RV data” worldwide. It is exactly the three RV series collections listed by the NASA Exoplanet Archive bulk-download page at the stated retrieval date. Additional observations remain with instrument teams, publications, and other archives. The included sample is not complete for stars without known planets and is not suitable for estimating intrinsic planet occurrence without a survey selection model.

Reported time coordinates include JD, HJD, BJD, MJD, UTC-labeled variants, and the MARVELS FCJD convention. The initial per-curve fits retain the reported time convention and do not combine curves across time standards. RV zero points, instrument configuration changes, stellar activity, secular acceleration, trends, and non-Keplerian signals may still produce model mismatch. The current model includes one jitter term per curve; this is not a full correlated stellar-noise model. The model-count prior and orbital prior bounds also affect $P(N\mid d)$ and will be varied in sensitivity tests.

Most importantly, the full corpus has not yet been sampled. No statement about a systematic change to the exoplanet census, multiplicity distribution, or individual published planet interpretation is supported until the full per-curve runs, reference matching, prior sensitivity, convergence assessment, and duplicate analysis are complete. The current deliverable is a data-audited analysis draft and executable preparation pipeline.

## 6. Data and code availability

PCAT source code is available at <https://github.com/tansudaylan/pcat>. The reproducible RV collection preparation script is `prepare_archive.py` in this directory. Raw archive files are retrieved into a user-selected external data directory and are not redistributed in this repository. The repository's current worktree includes unrelated uncommitted changes; a release DOI and exact source-tree archive must be created for the completed analysis.

## References

- Daylan, T., Portillo, S. K. N., & Finkbeiner, D. P. 2017, *ApJ*, 839, 4, [doi:10.3847/1538-4357/aa679e](https://doi.org/10.3847/1538-4357/aa679e).
- Howard, A. W., & Fulton, B. J. 2016, *PASP*, 128, 114401, [ADS](https://ui.adsabs.harvard.edu/abs/2016PASP..128k4401H/abstract).
- Thomas, D. et al. 2016, *PASP*, 128, 045003, [ADS](https://ui.adsabs.harvard.edu/abs/2016PASP..128d5003T/abstract).
- NASA Exoplanet Archive, contributed RV series and bulk-download snapshot, retrieved 2026-10-11, [archive documentation](https://exoplanetarchive.ipac.caltech.edu/docs/contributed_data.html).
