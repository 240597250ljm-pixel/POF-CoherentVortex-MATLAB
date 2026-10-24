# Changelog

## 1.0-rc3: method-level scientific correction

- Replace Liutex-intensity-weighted vortex centroids with geometric component centroids (matching final Fig05).
- Use original 4σ truncation for Liutex R field smoothing; retain 3σ for Gaussian coarse-grained flux and valid support.
- Record signed circulation and its magnitude separately.
- Make multi-frame spatial null-model RNG seeding safe: no default reseed per call; reject incomplete null ensembles rather than silently averaging failed placements.
- Add analytic identity tests for 2C horizontal flux (Galilean and cubic scaling), Liutex expansion invariance, S3 parity, core centroids and discrete Gaussian support.
- No experimental data or figure outputs included or claimed.

## 1.0-rc2

Original reference-method implementation.

### Public repository metadata (2026-10-24)
- Updated repository visibility wording and clarified that a public source snapshot does not yet grant an open-source reuse license.
- Scientific algorithms remain unchanged from v1.0-rc3 method-audited snapshot.
