# Lateral sensitivity and intrinsic spatial frequency data

Human psychophysical measurements accompanying the lateral sensitivity and intrinsic spatial frequency (ISF) analyses in Alexander A. Minns's doctoral thesis, *Wave-Based Cortical Coupling as a Mechanism for Spatial Frequency Selectivity and Gain Control* (Loughborough University, May 2026).

**Start with [the measurements guide](docs/MEASUREMENTS.md): raw staircases → reversals → thresholds → lateral sensitivity → ISF with error bars.** There is also a [short PDF guide](docs/Lateral_data_guide.pdf). For scientific reuse, read [methods](docs/METHODS.md), [known limitations](docs/KNOWN_LIMITATIONS.md), and the [data dictionary](data_dictionary.csv).

The main browsing tables are [staircase thresholds](data/processed/staircase_thresholds.csv), [lateral sensitivity](data/processed/lateral_profiles_recomputed.csv), and [ISF endpoints](data/processed/isf_endpoints.csv). The [endpoint overview](figures/measurements/isf_endpoints_with_percentile_ranges.png) shows the estimates with explicitly defined error bars. Fitting-code comparisons are secondary diagnostics.

This package contains **Subjects 1–4**. Subjects **1, 2 and 4** are the primary thesis ISF subjects. Subject **3** is included for completeness and is explicitly marked as diagnostic: the thesis excluded these ISF results from its main interpretation because of poor fit quality. Public IDs P01–P04 equal thesis Subjects 1–4.

## What is included

| Material | Contents |
|---|---:|
| Selected luminance conditions | 17 |
| Preserved numeric trial files | 560 |
| Recorded trials | 115,445 |
| Reconstructed staircases | 2,469 |
| Deterministic lateral profile points | 307 |
| Archived paired profile/spread inputs | 850,000 |
| Archived fit rows | 850,000 |
| Fits passing the saved RA FMS criterion | 280,785 |

The archived resampling inputs contain 15,350,000 spatial points. They are consolidated into compressed CSV files, with hashes that allow all 1,700,000 original two-column input files to be reconstructed and checked. Original trial and fit files retain their numeric bytes inside gzip files.

## Choose a starting point

- **See the lateral measurements:** `data/processed/lateral_profiles_recomputed.csv` and `figures/diagnostics/*_lateral_profiles.png`.
- **See the saved ISF results:** `data/processed/*_saved_fits.csv.gz`, `fit_retention.csv`, and `isf_mode_summaries.csv`.
- **See the supplied thesis figures:** `figures/reference/`. These four PNGs are unchanged.
- **Inspect trials:** `data/raw/trial_files.csv` identifies the preserved files in `data/raw/source_trials/`.
- **Inspect resamples:** `data/bootstrap/`, with one compressed table per condition.
- **Check or rerun analysis:** [reproduction instructions](docs/REPRODUCIBILITY.md).

## Verification and interpretation

All 17 deterministic profiles and their spreads reproduce the saved input-zero files exactly. All 15,349,693 nonzero-input spatial points are compatible with the same five-of-nine paired reversal-subset calculation within a tolerance of 1e-12. All 17 accepted-fit counts match thesis Tables 4.1 and 5.1, and all 13 primary-subject candidate-cluster count pairs match Table 5.2.

**Historical fitting is not fully reproduced.** Only 12 of 85 fixed-sample refits match both saved frequency and FMS to four decimal places with the supplied current helpers. The mapping from each historical fitted row to its original input is assumed for example replays, not established. Preserved results and newly computed diagnostic fits are clearly separated. See [known limitations](docs/KNOWN_LIMITATIONS.md) before interpreting a diagnostic as an original thesis result.

## Scope, citation and rights

This is a separate thesis lateral-sensitivity/ISF package. It does not duplicate the CSF/PSF data release, and it has no assigned release DOI or public repository URL. **Do not cite the CSF/PSF Zenodo DOI as the DOI for this package.** The `2024` directory label follows the requested project name; selected measurements include 2024 and 2025 sessions.

The full thesis and original scripts with local paths are not distributed here. Exact source hashes and extracted author-function provenance are included. Participant names, original source directory identities and acquisition dates are omitted from public metadata. This is a technical preparation, not an independent determination of consent for public sharing.

Licensing and contributor attribution need the owner's choice before public release; see `LICENSE_DATA.md`, `LICENSE_CODE.md`, and [release notes](docs/RELEASE_NOTES.md). No public release or upload has been made as part of preparation.
