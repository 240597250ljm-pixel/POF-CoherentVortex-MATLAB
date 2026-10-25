# Method-to-code alignment, October 24 manuscript/SI

This document reports *what the code computes*. `config/pof_paper_parameters.m` centralizes manuscript numbers; `pof.analyzeFrame` composes the frame-level reference methods. It does not assert scientific regression to a final production image or CSV. Required evidence for exact figure reproduction: frozen frame keys and raw velocity datasets, final producer script and parameters, frozen output CSV, and value-level comparison.

| Manuscript / SI method | Callable method | Alignment / caveat |
|---|---|---|
| PIV units, 0604 m/s; 0605 pixel + 50μs | `io/pof_read_vec_field`, `io/pof_parse_vec_header`, `io/pof_normalize_vec_units` | Header-based fail-closed conversion; 0605 114.910004 μm/pixel confirmed by experimenter but absent from header |
| Planar Liutex, 2×2 derivatives, 0.8-grid Gaussian R smoothing truncated at 4σ | `pof.planarLiutex` | Analytically tested solid rotation/shear; does not reconstruct 3D |
| P80 positive R, min 4 pixels in **vector grid**, 2-cell margin | `pof.detectVortices` | Discrete 8-connected components, hole fill, **geometric** equal-area centroid (as in the final Fig05 `regionprops`) |
| Γ = ∫ ωz dA and req = sqrt(A/pi) | `pof.detectVortices` | Field units m/s, coordinate spacing meters |
| Finite-window Voronoi and compact / large fractions | `pof.voronoiBounded` | All-pair half-plane clipping; check partition error; frame CV sample std |
| CSR_FULL / CSR_DETECT / SSI_CORE | `pof.voronoiNulls`, reference `finalize_POF_IIIB_Voronoi_from_Fig05_v1.m` | Caller should seed RNG ONCE per frame series; implemented with 12/24/24/12 replicates and 0.5/0.75 radius exclusion; the paper\'s random draws are not preserved, so its exact final null table is not regenerated |
| Gaussian coarse graining, stress-strain contraction | `pof.horizontalFlux` | 2C only; mean removal, default 0.8-grid U/V pre-smoothing (4σ kernel), separate ℓ-FWHM symmetric 3σ coarse graining and valid-domain crop |
| Structure function S3L, 14 bins | `pof.structureFunctionS3` | Unique half-plane offsets and pair-weighted bin means on one frame; frame medians/MBB are separate |
| Near field ρ≤1.5, BG ρ>2 | `pof.nearFieldMasks`, `pof.nearFieldStats` | Transition excluded, overlapping coverage union |
| R90 / A90 / fraction of negative | `pof.strongLiutexStats` | Uses positive R P90; measured area fraction not set to 0.1 |
| Scalar gradient/strain matching | `pof.matchGradientStrain` | Deterministic 6×6 quantile strata and standardized-log caliper, not causal |
| Frame serial dependence | `pof.movingBlockCI` | Circular moving-block of frames; 10/3000 default; no spatial-pixel IID bootstrap |

## Velocity pre-smoothing alignment with revised manuscript

Manuscript Sec. II E and SI Sec. 1.1 explicitly specify Gaussian velocity pre-smoothing (U,V sigma=0.8 vector-grid cells) after subtracting frame-wise mean, before separate ℓ-FWHM Gaussian coarse graining. The callable kernel **defaults to 0.8** for this reason; `opts.preSmoothSigmaCells=0` explicitly disables it for a sensitivity check. Pre-smoothing uses a discretely normalized symmetric ±4σ kernel; scale-dependent coarse graining uses ±3σ support and 3σ valid-domain exclusion. Gaussian smoothing of the Liutex R field is a separate operation. This method alignment does not certify numerical replication of final figures or resolve all historical producer-data provenance.

## Data and external scripts

This archive contains no 0604/0605 raw VEC/DAT, accepted-frame CSVs, final flux CSVs, Voronoi null replicates, or manuscript figure PNG/TIFF/PDF. A fixed 81773 vortex catalog and a separate 81775 historical extraction occur in the research history. The formal manuscript primary catalog is 81773 and 2935 retained frames; do not merge them. Old scripts are separately archived. Figure mapping in `docs/FIGURE_PROVENANCE.csv` states whether each output can be exactly recreated from this distribution.

**Quantile convention:** this release uses explicitly programmed inclusive linear interpolation `r=1+(n−1)p/100`, matching the original `finite_percentile` method, instead of MATLAB-version-dependent `prctile` defaults.

## Static scientific audit correction (v1.0-rc3)

- `detectVortices`: centroid is geometric (unweighted) as in the final Fig05 production `regionprops` call; previous release wrongly weighted with Liutex intensity. `Gamma_m2ps` is signed; new `AbsGamma_m2ps` is the main manuscript population metric.
- `planarLiutex`: R smoothing at sigma=0.8 cells is truncated at ±4σ, as in final Fig04/05 source; `horizontalFlux` coarse-graining still truncates at ±3σ and excludes the corresponding edge strip.
- `voronoiNulls`: default uses the ambient RNG state, rather than fixing the same seed for every frame; seed once per six-condition experiment, not once for each individual frame.
- `horizontalFlux`: revised manuscript and SI explicitly state 0.8-grid U/V pre-smoothing, now enabled by default; disabling it remains an explicit sensitivity option.
