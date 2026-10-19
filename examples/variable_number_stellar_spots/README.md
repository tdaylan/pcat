# Variable-number stellar spots

How many spots can a light curve support? This example simulates 25.6 days of
Poisson photometry from a star with three rotating dark spots, then asks PCAT to
fit catalogs containing zero to five spots. The rotation period is fixed at 3.2
days. Each spot has a sampled depth, phase, and width. The frame sequence shows
the observed counts, a PCAT model sample, and their residuals.

Run a short pipeline check with:

```bash
python examples/variable_number_stellar_spots/variable_number_stellar_spots.py --smoke
```

The simulated data are illustrative and are not observations of a specific
star.