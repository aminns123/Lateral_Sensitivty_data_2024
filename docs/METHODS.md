# Methods and source interpretation

## Evidence and scope

The preparation used the supplied thesis, storage driver, current fitting driver, current plotting driver, both author helper libraries, selected raw data, saved resampling inputs and saved fit files. See `config/source_inventory.json` for source hashes and `config/function_provenance.json` for exact extracted-function hashes and original line numbers. This archive covers empirical lateral sensitivity and ISF, not all theory, simulations or alternative model comparisons in the thesis appendices.

Relevant thesis locations use printed page numbers, followed by PDF page numbers: experiment fit exclusion Table 4.1, p.123 (161); flanker settings Table 4.2, p.126 (164); reversal subsets Eqs.4.14–4.16, p.136 (174); main ISF screening pp.150–152 (188–190); accepted counts Table 5.1, p.156 (194); cluster counts Table 5.2, p.161 (199); models Eqs.B.19–B.25, pp.220–221 (258–259); reversal resampling pp.223–225 (261–263); spread and likelihood pp.228–230 (266–268); outlier methods and FMS pp.231–233 (269–271). Appendix A provides the model background rather than extra empirical observations for this package.

## From trials to lateral sensitivity

The historical loader groups rows within a source file by absolute probe position and staircase ID, retains the first 200 rows in a group and uses the signed position from its first row. The reconstructed table records the number of distinct signs to make this behavior inspectable. Original row indices are one-based; analysis and reversal indices are zero-based.

The source digital contrast calculation is `(probe_intensity - background_intensity) / (1 - background_intensity)`. These are normalized digital-intensity quantities; the preparation does not retrospectively convert trial values using an unverified physical calibration. Reversals are detected using the exact author helper `count_reversals_HighLow`.

For deterministic input zero, each staircase threshold is the geometric mean of its final five detected reversal contrasts. Thresholds are combined by their arithmetic mean across included staircases at a position. Sensitivity is the reciprocal of that mean threshold. The analysis quantity is `R(x) = ln(S_flanker(x) / S_base(x))`.

Spread is the sample standard deviation, with denominator n-1, of all available pairwise staircase-level `ln(T_base / T_flanker)` values at the position. It is not a standard error. Supplied plotting conventions use half that spread as the error-bar half-width. The weighted fit uses the full spread with a floor of 0.05. Neither is a population confidence interval.

The source directory enumeration index is preserved because recorded staircase exclusions depend on it. All trials, including excluded staircases, remain in the archive. `validation/exclusion_candidate_checks.csv` compares no exclusions, recorded exclusions and the active Subject-1 storage-script exclusions. `config/staircase_groups.json` marks the candidate used to match input zero. This inference is supported by exact profile/spread agreement and by compatibility of every archived nonzero profile point with its reversal-subset pool; it does not establish why each historical exclusion was made.

## Reversal-subset resampling

Five distinct indices are selected from reversal indices 1–9 (human-numbered reversals 2–10), giving 126 possible subsets. At a probe position, the same subset indices are applied to each retained staircase and to both base and flanker conditions. Thresholds are computed within each staircase before averaging across staircases. Base and flanker trials came from separate sessions; pairing here is computational.

Every selected condition has 50,000 archived input profiles and 50,000 archived fit rows. Input index 0 is the deterministic last-five profile; indices 1–49,999 are resampled profiles. The current storage routine generates unique joint vectors of subset choices across positions; its historical random state was not saved. The new seeded utility uses explicit modern random state and unique joint vectors, and writes only outside this archive. It does not claim to recover historical ordering or add a historical fit-row mapping.

`audit_resampling.py` checks every nonzero-input spatial point against the complete set of 126 possible jointly calculated ratio/spread pairs for its position, using the baseline-matching exclusion candidate. All 15,349,693 points are compatible at absolute tolerance 1e-12. This establishes pointwise computational compatibility, not uniqueness of the generating subset or the historical random sequence.

## Stored fits and fit screening

The selected source branch is the five-reversal geometric-mean, all-conditions, two-residual, `MLE_greenFUNC_phase` analysis. Each saved file has 50,000 rows and 17 columns; numeric precision and tokens are retained. The first 16 column roles are interpreted from source/plot usage; the auxiliary seventeenth field is preserved with a qualified name. The historical writer itself is not version-matched to these files.

Residual-analysis (RA) FMS is original column 0; RA intrinsic frequency is column 9, both zero-based. Inclusion uses `0.8 <= FMS_RA <= 1.0`. No historical optimizer-success flag exists in the saved file; `fit_success_recorded` is blank, not True. LOF and RANSAC alternative results are preserved but do not drive the principal screen. FMS compares normalized amplitude spectra, not phase agreement. Its use on finite and potentially irregular position sequences is a diagnostic inherited from the author helper, not a guarantee of frequency identifiability.

The spatial model is thesis Eq.B.25: `A exp(-lambda |X|) [cos(k X + phase) - sign(X) sin(k X + phase)]`, with `X = x + 0.5` in degrees and `f_n = k / (2*pi)` in cycles per degree. It is distinct from the independently weighted sine/cosine model Eq.B.19.

## Candidate modes and diagnostic examples

The current peak helper locates extrema on a smoothed 31-bin histogram spline (smoothing parameter 0.02), although it also computes a KDE. The separate display/mass KDE uses Scott bandwidth multiplied by 1.0, evaluated at 1,000 points over the accepted frequency range. Candidate intervals are selected with the exact author helper. Any median fallback is recorded in `peak_locator`; no fallback is silently described as an original peak.

Two percentile variants are reported: 0.5–99.5 from the thesis/figure legend, and 2.5–97.5 from the current driver. Modes require integrated KDE mass at least 0.10 inside their percentile interval; at most two are displayed. Candidate IDs follow increasing frequency. A candidate can have a nonzero count and still fail the display rule. The largest retained mass defines the dominant mode independently at each luminance.

Representative diagnostics separately replay the source plotting index rule and the maximum saved FMS within each candidate cluster. They assume saved row n corresponds to input n-1. Fits use the supplied current helpers, ordered parameter names, and the source filtered-residual indexing behavior, with outcomes and warnings recorded. The original four reference PNGs remain untouched. See the limitations document for differences from the prose description and historical saved fits.
