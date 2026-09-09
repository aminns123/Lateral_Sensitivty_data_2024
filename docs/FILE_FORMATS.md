# File formats and provenance

`data_dictionary.csv` documents every distributed CSV column and every column of the headerless original numeric formats. Booleans use True/False, indices explicitly state their origin, and blank numeric/status values mean unavailable rather than zero. JSON null also means unavailable. Decimal points are used throughout.

## Original files and de-identification

`source_manifest.csv` records each preserved standalone numeric file or supplied reference image. A gzip-wrapped source file must be decompressed before checking its original SHA-256. Raw numeric content is unchanged. Public IDs replace original directory/participant identities and session filenames. The captured source-directory index remains available for matching exclusions. Original metadata text and the private directory-to-public-ID map remain in a separate owner-only provenance archive, not in this repository.

`data/raw/stimulus_metadata.csv` contains retained source key/value strings, not a claim that every parameter snapshot described the executed acquisition. Names, dates and machine paths were removed from this public metadata. Unrecognized raw columns remain explicitly uninterpreted; they were not dropped or recoded. Source column positions are zero-based in the dictionary.

## Consolidated bootstrap inputs

Each `data/bootstrap/*_profiles.csv.gz` contains input indices 0–49,999 and one row per spatial point. Ratio and spread share the recorded position token. `validation/*_bootstrap_sources.csv.gz` records original byte counts, hashes, row counts, line-ending style and final-newline presence for each of the two original files per input. These records let the verifier reconstruct the exact original two-column tab-separated text without distributing 1.7 million tiny files.

`data/source/baseline_profiles` separately preserves input zero. `data/source/example_inputs` preserves the 139 paired inputs used in the fixed sample and representative replays. These are intentional duplicates for convenient small-example use, not extra independent data.

## Configuration JSON

- `conditions.json`: one record per condition, public subject IDs, original label, source drive settings, precise thesis table values and main-analysis inclusion. Null thesis table values mean the subject was not listed in that table.
- `analysis_settings.json`: explicit screening, resampling, coordinate, uncertainty and peak-summary conventions.
- `staircase_groups.json`: condition, experiment, captured file index, public source ID, original within-file staircase ID, public staircase ID, signed position, contrast reversals, one-based source rows and the baseline-candidate inclusion flag. This is generated from preserved trials.
- `representative_requests.json`: condition, candidate ID, selection rule, assumed input index and spline peak used to initialize a diagnostic replay.
- `source_inventory.json`: source script/thesis identities by basename and SHA-256, without local paths. The two helper files share a basename but have distinct source IDs/hashes.
- `function_provenance.json`: distributed helper module, function name, source ID, original start/end lines and exact extracted function-text hash. Imports were adapted; function bodies were preserved.

## Checks and missing information

Validation JSON files summarize executable checks; accompanying CSVs give condition- or example-level evidence. `fit_success_recorded` is deliberately empty for every historical fit because the source file did not save it. `optimizer_success` applies only to a new replay. `assumed_input_index` is deliberately distinguished from the actual index recorded in archived input files.

The final `manifest-sha256.csv` authenticates the distributed file bytes and excludes itself. `source_manifest.csv` authenticates decompressed historical source bytes. These serve different purposes.
