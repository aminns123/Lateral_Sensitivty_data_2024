# Reproducing and checking results

The public scripts use paths relative to their own location. They do not require the preparer's computer paths. Keep original reference files intact; perform reruns in a working copy. Commands below run from the repository root using Python 3.12.

## Install dependencies

Create an isolated Python environment, activate it, then run:

```text
python -m pip install -r requirements.txt
```

Versions in `requirements.txt` are the preparation environment, not a recovered historical environment. `verify_package.py` uses NumPy and standard-library modules; the refit/mode/plot scripts additionally use SciPy and Matplotlib. No scikit-learn dependency is needed for the selected RA replay; archived LOF/RANSAC results are preserved without re-estimating them.

## Verify the downloaded reference package

```text
python code/verify_package.py
```

This checks the complete manifest, source hashes, trial dimensions, all 850,000 fit-row conversions and flags, thesis count comparisons, exact deterministic profile/spread agreement, extracted function bodies and dictionary coverage. It reconstructs all 1,700,000 original profile/spread files in memory from the compressed tables and verifies their byte counts and hashes. It does not write over data. Allow several minutes for the complete check.

`validation/integrity_report.json` records the preparation check made before the final self-excluding manifest was written. Its manifest-file count is therefore zero; the final manifest was subsequently checked separately. The SHA-256 manifest excludes only itself, `.git`, `.venv` and `__pycache__` files. Unlisted files cause verification to fail. Gzip headers are not historical evidence: decompressed source bytes are what `source_manifest.csv` authenticates.

## Recompute analysis in a working copy

```text
python code/reproduce_raw.py
python code/summarize_modes.py
python code/refit_examples.py
python code/plot_data.py
python code/audit_resampling.py
python code/export_measurements.py
```

These overwrite generated tables/diagnostics in that working copy, never raw/source data. Recomputed gzip headers or library-dependent numerical outputs may change the original full-package manifest. After rerunning, use `python code/verify_package.py --skip-manifest` to check the scientific invariants and source preservation. Preserve the original reference download for byte-level comparison.

`reproduce_raw.py` extracts staircases/reversals and compares explicit exclusion candidates to the saved deterministic input zero. `summarize_modes.py` reads the preserved fits and produces both percentile variants. `refit_examples.py` reruns the fixed audit sample and both representative-selection rules, recording success, errors, parameters and removed points. It is a current-code diagnostic, not a complete regeneration of the historical fit files. `plot_data.py` reads these tables to draw labelled figures. `audit_resampling.py` checks every nonzero input point against its 126 possible paired reversal subsets.

`export_measurements.py` provides the primary browsing path: used trial rows with reversal flags, staircase thresholds with inclusion flags, a compact endpoint table, staircase examples and the endpoint overview. Run it after the raw and mode summaries above. It uses saved fits for the endpoints and does not perform new fitting.

## Make a new, seeded resampling run

Use an output directory that does not already exist and is outside the reference repository:

```text
python code/resample_profiles.py --condition P01_L010 --seed 20260908 --samples 1000 --output ../lateral_new_resamples
```

The output records new sample IDs, chosen subset indices and a seed. It uses the baseline-matching exclusions and shared reversal indices at each position. It creates new profiles only, with no deterministic input-zero entry and no historical fit mapping. No historical data are replaced. The package intentionally does not supply a purported historical full-refit command when the writer/version identity is unresolved.

## Read large tables in Python

```python
import csv, gzip
with gzip.open('data/bootstrap/p01_l010_profiles.csv.gz', 'rt') as f:
    for row in csv.DictReader(f):
        # Each row is one position in one preserved input profile.
        if row['source_input_index'] == '0':
            print(row)
        else:
            break
```

CSV is UTF-8 with comma separators and decimal points. `.tsv.gz` files are gzip-wrapped original numeric files; read them with whitespace splitting or `numpy.loadtxt`. They have no headers. See `docs/FILE_FORMATS.md` and `data_dictionary.csv`.
