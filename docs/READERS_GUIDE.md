# A guide to this folder

You can explore the results without running Python. Begin with `figures/reference`, which contains your four supplied images, and `figures/diagnostics`, which contains clearly labelled new views of the data. The ordinary `.csv` summary tables open in a spreadsheet application. A `.csv.gz` file is a compressed CSV: extract it first or read it directly with Python. The largest tables may exceed a spreadsheet application's worksheet capacity, so use the smaller summaries for browsing.

## The measurements in plain language

The participant judged whether a small probe appeared above or below fixation. Repeated trials changed the probe intensity to estimate the contrast needed for detection. Measurements were made both with and without an indirect flanker.

**Lateral sensitivity** describes how detection sensitivity changes with probe position. Sensitivity is the reciprocal of the contrast threshold. The stored log ratio compares sensitivity with the flanker to sensitivity in the base condition: positive means greater sensitivity with the flanker, negative means lower sensitivity, and zero means equal sensitivity.

**Intrinsic spatial frequency** is inferred by fitting an oscillating, decaying curve to this log ratio. It is a model-derived spatial scale, not a direct measurement of a cortical wave. Repeating the calculation with different reversal subsets produces many possible fitted frequencies. Their distribution describes the stability of this estimator for these measurements.

## How the folders connect

1. `data/raw` contains the preserved trial records and de-identified stimulus settings.
2. `data/processed/staircase_reversals.csv` identifies the reversals detected in each staircase.
3. `lateral_profiles_recomputed.csv` combines the last five reversals to reproduce the original deterministic profiles.
4. `data/bootstrap` contains the original resampled profiles and their spreads. There are 50,000 inputs per condition, including deterministic input zero.
5. `data/source/saved_fits` preserves the original fitted results; the corresponding files in `data/processed` add column names and screening flags without changing the saved numeric tokens.
6. `isf_mode_summaries.csv` describes candidate frequency clusters, their percentile ranges, and whether they pass the display rule.
7. `validation` records the checks. `code` lets a reader recompute profiles, summaries, diagnostic fits and figures, or make new seeded resamples.

The filename `p01_l010` means thesis Subject 1 at the archived luminance label 10 cd/m². Labels remain as originally archived. More precise thesis values, and a disagreement at Subject 1's label 20, are recorded separately in `config/conditions.json`.

## Which figures should I use?

The four images in `figures/reference` are your originals. New `*_lateral_profiles.png` images show the independently reconstructed last-five-reversal measurements. The other diagnostic figures replay two different rules for choosing representative examples. They are labelled as new fits and should not silently replace the submitted or thesis figures.

Both representative rules rely on an assumed relationship between fitted row numbers and input numbers. The current plotting code's indexing differs from selecting the highest saved FMS in a cluster. We retained both views so that the difference can be examined without losing the original figures.

## How much has been checked?

The original lateral profiles reproduce exactly. The entire archived resampling collection is compatible with the documented reversal-subset method. Saved acceptance counts and primary-subject cluster counts agree with the thesis. Original bytes are checked with cryptographic hashes.

There is a remaining historical fitting issue: the current scripts do not recreate most of the sampled historical fitted results. The stored fits are preserved; the guide does not present a refit as if it were the original. Full details are in `KNOWN_LIMITATIONS.md`.

## What remains before publishing?

The local package is ready to inspect and use. Before a public release, choose a data licence and code licence, confirm the intended contributors and public-sharing permissions, and resolve or explicitly retain the documented scientific qualifications. Then create a new GitHub repository and, if desired, its own Zenodo record. The DOI for the CSF/PSF repository belongs to that other dataset.
