# Known limitations and decisions

These qualifications distinguish source evidence, reproducible results and unresolved historical details. They are part of the dataset, not hidden preparation notes.

## Historical fitting and representative figures

- The supplied current MLE-writing script does not match the selected historical 17-column files. Its current branches include incompatible output-width and parameter-name constructions. It was inspected, not executed against the original data.
- A fixed audit refitted inputs 0, 1, 1234, 25000 and 49999 in all 17 conditions with the current helpers: **85 successful optimizations, zero execution errors, but only 12/85** matched both saved RA frequency and RA FMS to four decimal places. A successful optimizer is not proof of historical reproduction.
- Equal input/fit counts do not establish their original row-to-input mapping. Example fits use row n to input n-1 as an explicit assumption. The archived inputs and fit rows remain independently identifiable.
- The supplied representative selection function indexes a filtered score vector using original row IDs, then uses a fallback that can skip candidates. It differs from the highest saved FMS within a cluster in **all 27 retained-mode selections** audited here, including Subject 3 diagnostics. Both rules were replayed (54 successful refits); neither set is asserted to be pixel-identical to the four reference figures.
- The residual helper returns a filtered residual array. The supplied caller then applies indices of the two largest values in that filtered array to the original profile. Diagnostic replays preserve this indexing behavior and record the removed indices. It may not remove the intended two original observations. A corrected scientific analysis would need separate outputs and renewed interpretation, rather than overwriting saved results.
- The portable refit wrapper supplies parameter names in a stable order; the original fitting script uses a set. Only needed helper functions are extracted, with import changes to avoid namespace collisions. Their function bodies are hash-verified.

## Thesis versus current plotting settings

- The thesis describes locating peaks on a KDE; the supplied peak helper actually locates them on a smoothed histogram spline. The current helper nevertheless exactly reproduces all 13 primary-subject cluster-count pairs in Table 5.2.
- The thesis/figure legend uses 0.5th–99.5th percentiles, while the current driver uses 2.5th–97.5th. Both summaries are included. New diagnostic panels explicitly use the thesis percentile endpoints with the supplied histogram-spline peak locator.
- Reference figure captions describe median markers and highest-FMS profiles. New diagnostics label their spline peaks and selection rule explicitly. They should not be substituted for original figures without deciding which scientific convention is intended.
- Table 5.2 describes retained clusters, but some listed candidate clusters fail the stated KDE-mass display threshold. Counts and `retained_for_display` are separate fields here. A count of one does not constitute evidence for a retained mode.
- Appendix B contains illustrative simulations and examples beyond the empirical release. Their plotted counts need not equal the final Chapter 5 selection. This archive does not claim to reproduce those simulations or every appendix figure.

## Acquisition metadata and labels

- Subject 1's archived luminance label 20 is listed as **15.0 cd/m² in Table 4.2**, but 20 in the results and supplied figures. Neither value is silently corrected. Separate fields preserve the archived label and thesis table values.
- Subject 2's label 26 corresponds to thesis 26.15; Subject 4's label 10 corresponds to thesis 10.2. Subject 3 is absent from Table 4.2; those table-specific fields are null.
- Older Subject-1 flanker folder settings such as 100 or 150.3 use a legacy scale; they are not cycles per degree. Thesis Table 4.2 gives the physical frequencies. A single conversion must not be applied to all conditions, because later folder values already use the smaller scale.
- Some Subject-1 source flanker contrast settings differ from the thesis's reported 0.018. Digital-drive settings and measured display contrast need not be identical. The precise calibration relationship was not reconstructed, so original settings and thesis statements remain distinguishable.
- Stored parameter snapshots contain differences from the thesis, including stimulus duration 250 versus 300 ms and staircase step settings 0.3/0.25245 versus 0.36/0.30294. Some metadata may be defaults rather than a record of executed trials. The archive preserves them and does not treat them as proof of actual timing.
- Recorded staircases are not uniformly ten-reversal sequences: 2,030 have 10, 429 have 11, eight have 12, and two have nine detected reversals under the author helper. The two nine-reversal staircases are excluded by the matched baseline candidate. The archived point count also varies (12–28 per condition), rather than uniformly following the 16/20-position description.
- Exclusion indices originally depended on local directory order, including subdirectories. The captured indices are preserved. Current Subject-1 code applies a global exclusion set across luminances; its scientific rationale cannot be recovered merely from matching numbers.

## Limits of inference and release status

The reversal-subset distributions describe within-observer estimator stability. They are not distributions of independent participants or population confidence intervals. Mode structure alone does not establish distinct neural states. Subject 3 remains excluded from the primary thesis ISF interpretation despite inclusion of its data here.

Historical seeds, original optimizer diagnostics, a version-matched fitting environment, complete calibration records, and a certified historical row-to-input linkage were not recovered. The archive makes strong claims about preserved bytes, raw-profile reproduction, resampling compatibility and reported counts; it makes no claim that every historical fit has been regenerated.

The rights holder still needs to choose publication licences, confirm contributors and ensure public sharing is consistent with the study's permissions. No DOI from another dataset applies to this package.
