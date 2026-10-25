 Planar vortex organization and horizontal interscale transfer (MATLAB)

**Version 1.0-rc4 — manuscript-aligned method reference / public source snapshot.**

This MATLAB package provides a **unit-safe, modular reference implementation** of the principal methods described in the October 24, 2026 revision of *“Coherent-vortex reorganization and near-field localization of horizontal interscale-transfer activity in buoyancy-driven free-surface turbulence”* (*Physics of Fluids* manuscript). It intentionally does **not** claim that the generated code has reproduced the frozen paper figures. This is a methods audit, not a guarantee that every experimental or historical implementation is free of error. Scientific numerical regression using the authors' raw data remains separate from synthetic method tests.

## Quick start

Use MATLAB R2022b+ with Image Processing Toolbox. In the repository root:

```matlab
run_pof('status')    % informational
run_pof('methods')   % supported reference functions
run_pof('test')      % mathematical synthetic tests
run_pof('example')   % synthetic 57 x 57 field, no raw input
```

Optional, **read-only** data-path check: `run_pof('preflight')` after setting environment variables described in `config/pof_config_local.m`. Preflight is **not** required to upload the repository.

## Actual supported public methods (`+pof/`)

| Function | Scientific definition |
|---|---|
| `analyzeFrame` | Safe one-frame method pipeline from calibrated 57×57 physical U/V grids, not full raw-data QC |
| `planarLiutex` | Nonnegative 2×2 planar diagnostic, not the 3D Liutex vector; Gaussian smooth R at sigma=0.8 grid cells with the original 4σ numerical kernel cutoff |
| `detectVortices` | Positive-R P80; hole fill; 8-connected areas ≥4 cells; exclude components intersecting 2-cell frame margin; geometric centroids; Γ, req |
| `voronoiBounded`, `voronoiNulls` | Bounded rectangular Voronoi; CSR_FULL, CSR_DETECT and SSI_CORE 0.5/0.75 conditional null references; initialize the RNG once per study |
| `horizontalFlux` | Two-component Gaussian coarse-grained Πℓ = −τij Sij; symmetric extension; exclude 3σ edge strip; **not complete 3D energy flux** |
| `structureFunctionS3` | Half-plane unique integer offsets; 2-grid minimum, 1/4-field maximum; 14 radial bins; frame-level S3L |
| `nearFieldMasks`, `nearFieldStats` | ρ≤1.5 vs ρ>2 (transition excluded); activity; negative-event contribution/probability/conditional magnitude identity |
| `strongLiutexStats` | P90 of **positive** planar Liutex over valid region; true A90, F90 and enrichment (SI S4) |
| `matchGradientStrain` | Within-frame matching by log gradient/strain magnitudes; reports coverage/balance, not a causal estimate |
| `movingBlockCI` | Circular 10-frame/3000-draw moving-block bootstrap of *frame* statistics, not individual pixels |

Units: 0604 VEC U/V are already m/s; 0605 U/V are pixel displacements per 50μs, using the experimenter-confirmed 114.910004μm/pixel calibration (2.29820008 m/s per pixel); the original 0605 VEC header lacks spatial calibration. `io/` blocks implicit conversions.

### Deliberately explicit preprocessing

Manuscript Sec. II E and Supplementary Methods Sec. 1.1 specify two distinct steps after subtracting the frame-wise spatial mean: **velocity pre-smoothing** of U,V at sigma=0.8 vector-grid cells (discrete normalized symmetric kernel, ±4σ support), followed by **scale-dependent Gaussian coarse graining** at FWHM ℓ (±3σ support with 3σ valid-domain exclusion). Both `pof.horizontalFlux` and `pof.analyzeFrame` now default to the stated 0.8-cell pre-smoothing. Specify `opts.preSmoothSigmaCells=0` only to study the effect of omitting that preprocessing. Separately, Liutex vortex identification smooths the computed R field at sigma=0.8 cells; that operation is independent of velocity pre-smoothing. This is a documented **method-level reference**, not a certified numerical reproduction of all paper figures. See `docs/METHOD_ALIGNMENT.md`.

## Structure and provenance

- `+pof/`: runnable, unit-independent scientific reference kernels.
- `io/`: native Tecplot VEC reader, tested metadata interpretation, strict physical unit conversion; read-only local preflight.
- `config/`: local data-root configuration and fixed manuscript method parameters; raw PIV files not distributed.
- `reference/`: three **existing source scripts** for figures / fixed-catalog analysis; they are *not* run by `run_pof`, may use old fixed paths, and are marked provenance references only.
- `tests/`: synthetic mathematics checks, no experimental results.
- `docs/`: precise methodology alignment, missing data, figure claim boundaries.

Other 47 research/history source scripts from v1.0-rc1 are retained in a **separate legacy archive**, not silently converted into published production code. Historical code includes known m/s double conversions, independent-frame bootstrap and previous catalog versions. It is not part of the current callable API.

## Scientific and release boundaries

- Nominal Ra_L values: [0.463 0.926 1.39 1.85 2.32 2.78]×10^9; paper primary cohort: 2,935 frames (498/499/499/492/470/477), 81,773 detected vortices.
- Paper Gaussian FWHM ℓ: 11.03, 15.01, 20.43, 27.80, 37.83, 51.48 mm. Baseline near field: ρ≤1.5; 1.5<ρ≤2 excluded; background ρ>2.
- Inverse-cascade claims, 3D energy-flux claims, causality of vortices, full 3D vortex topology and absolute small-displacement PIV precision are **not** supported by this package.
- The software license and coauthor permission are **not yet decided**; see `LICENSE_PENDING.md`. The repository is currently public; a software reuse license and any further permissions remain pending (see `LICENSE_PENDING.md`).

**This public GitHub source repository provides inspectable reference routines, not certified end-to-end reproduction of the manuscript figures.**

**Kernel convention:** final Fig04/05 Liutex R smoothing uses a normalized ±4σ finite kernel; Gaussian coarse-graining for Πℓ uses ±3σ with a 3σ valid-domain exclusion. Vortex centers are geometric centroids, not intensity-weighted. See `docs/METHOD_AUDIT.md`.

**Quantile convention:** this release uses explicitly programmed inclusive linear interpolation `r=1+(n−1)p/100`, matching the original `finite_percentile` method, instead of MATLAB-version-dependent `prctile` defaults.
