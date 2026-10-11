# PCAT reanalysis of public radial velocities

This paper project targets the question: how does a transdimensional Keplerian-catalog posterior compare with the single catalog interpretation attached to each public RV curve?

## Data boundary

“All public RV data” is not a finite global archive: instrument teams publish additional series through papers, mission archives, and observatory repositories. The reproducible first corpus is therefore the complete set of three RV time-series collections listed by the NASA Exoplanet Archive bulk-download page as of 2026-10-08:

- `RADIAL`, the literature RV-curve manifest. The archive lists 1,072 download entries, including one repeated filename. The downloaded snapshot has 1,071 unique curves for 511 archive UIDs.
- `MARVELS_RV_bulk_download.tar.gz`, 6,320 curve files and 3,169 header `STAR_ID` values.
- `Fulton_rv_data.tar.bz2`, 76 stellar RV files with four documented instrument codes.

The collections are not a homogeneous survey or a statistically representative planet sample. The RADIAL files cover 10 time-reference labels and 68 instrument strings. Some curves for the same target and reported epoch conflict in velocity or uncertainty. MARVELS includes non-finite RV/error rows and repeated target/file segments. These data products will not be concatenated blindly or counted as independent system-level evidence. The preparation script preserves the original file, source reference, frame, instrument metadata, invalid rows, and repeated-epoch flags.

## Prepare the data

Use a directory outside the repository for the downloads. The script fetches the official RADIAL manifest and all unique files referenced by it, plus the MARVELS and Fulton bundles, extracts them safely, and writes normalized observations, per-curve metadata, and a quality summary:

```bash
python docs/papers/rv_catalog_reanalysis/prepare_archive.py \
  --data-dir /path/outside/the/repository/pcat-rv-archive
```

Raw files are not included in Git. Record the generated `audit_summary.json`, archive hashes, retrieval date, and any later curation decisions with each analysis.

## Run per-curve fits

```bash
python docs/papers/rv_catalog_reanalysis/run_curve_fits.py \
  --data-dir /path/outside/the/repository/pcat-rv-archive \
  --numbswep 100000 \
  --seed 17
```

Result tables and PCAT caches are isolated by seed, sweep count, planet cap, and minimum observation count, so a pilot cannot be resumed as though it were a differently configured run. Use a new seed for an independent chain.

## Analysis unit and interpretation

Each published curve file is the initial analysis unit. PCAT's RV likelihood fits Keplerian catalog size and orbital parameters, analytically marginalizes linear instrument offsets, and numerically marginalizes a jitter grid. All reported parameter posteriors remain per curve unless measurements have been shown to be non-overlapping, time-standard compatible, and attributable to the same stellar component. Cross-curve combinations require a separately justified duplicate and instrument-calibration model.

The eventual comparison will match each curve's source reference to the planet solution(s) reported in that reference and compare them with the PCAT posterior over zero to five Keplerian signals. Curves with too few valid epochs, unresolved same-epoch conflicts, or ambiguous target matches will be identified, not silently discarded. Aggregate planet-frequency claims require an explicit survey selection function and are outside the current inference target.

## Current status

The archive was audited and 7,051 of 7,467 curves meet the current runner's minimum-data and repeated-epoch eligibility rules. HD 142 run-length pilots remain unconverged: two 100,000-sweep chains give similar overall planet-count posteriors but show within-chain shifts, while 2,000-sweep runs vary by configuration. No revision to the published interpretation is supported. The complete per-curve inference sweep, archive-reference matching, sensitivity analysis, and final manuscript results are not yet complete.

## Sources

- NASA Exoplanet Archive, [bulk data download](https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/) and [contributed data sets](https://exoplanetarchive.ipac.caltech.edu/docs/contributed_data.html), retrieved 2026-10-11. The archive lists RADIAL as 1,072 published RV curves in 3 MB in the 2026-10-08 bulk listing.
- Thomas et al. (2016), MARVELS radial-velocity data, [ADS record](https://ui.adsabs.harvard.edu/abs/2016PASP..128d5003T/abstract).
- Howard and Fulton (2016), RV survey completeness and released observations, [ADS record](https://ui.adsabs.harvard.edu/abs/2016PASP..128k4401H/abstract).
- PCAT software, [GitHub repository](https://github.com/tansudaylan/pcat).
