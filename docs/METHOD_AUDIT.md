# Scientific method and definition audit (v1.0-rc3)

**Scope:** static, line-level comparison of callable MATLAB methods in `+pof/` against the October 24 revised manuscript and supplementary methods (both dated October 24, 2026), plus relevant final Fig04/Fig05 production-source snippets. Does not require proprietary raw PIV files and does not assert exact frame-by-frame production output.

| Method | Mathematical definition / correspondence | Status |
|---|---|---|
| Planar Liutex | `beta=(dv/dx-du/dy)/2`; `alpha=0.5*hypot(du/dx-dv/dy, du/dy+dv/dx)`; `R=2 max(|beta|-alpha,0)` only for `beta²>alpha²`; `lambda_ci=0.5 sqrt(max(4 det(G)-trace(G)^2,0))` | Checked, planar only; R smoothing ±4σ as in Fig04/05 producer |
| Vortex core | Gaussian R smoothing σ=0.8 grid cells, framewise P80 of positive R, 8-connectivity, hole filling, size ≥4 grid cells, 2-cell border exclusion | Checked; distinct from velocity smoothing |
| Vortex centroid, Γ, req | **Unweighted** geometric centroid of each retained core; signed Γ = Σ(vx−uy) ΔA and population |Γ|; req=√(A/π) | Corrected to final Fig05 producer; original rc2 wrongly intensity-weighted center |
| Rectangular Voronoi | Half-plane clipping against finite rectangular domain; `a_i=A_i/mean(A)`, frame CV sample stdev/mean, fractions `a<0.5` and `a>2` | Checked; all-pair clipping, no periodic wrap |
| Null models | CSR_FULL, CSR_DETECT and sequential inhibition `d_ij≥f(ri+rj)` at `f=0.5,0.75`; 12/24/24/12 realizations per frame | Checked definitions; default no repeated per-frame seed reset, minimum successful replicate gates (12/24/22/10); null-replicate averaging and framewise contrasts are downstream analysis, not claimed in one-frame kernel |
| Horizontal interscale flux | `tau_ij=overline(u_i u_j)−ū_i ū_j`; `S_ij=0.5(∂_j ū_i+∂_i ū_j)`; `Π=−(tau_xx S_xx+2 tau_xy S_xy+tau_yy S_yy)` | Checked 2C identity, symmetric 3σ convolution, 3σ edge exclusion; not full 3D cascade |
| Gaussian scale | `σ=ell/(2√(2 ln 2))`; ell in meters, **FWHM**; no hidden Hann window or solenoidal projection | Checked |
| Third-order S3L | `S3(r)=<[(u(x+r)−u(x))·rhat]^3>`; unique half-plane integer offsets, 2-cell min, quarter-domain max, 14 radial bins, pair-weighted means per frame | Checked code-level discrete convention |
| Near field and statistics | `ρ(x)=min_i |x−c_i|/r_i`; NF `ρ≤1.5`, transition `1.5<ρ≤2`, BG `ρ>2`; `A=E|Π|`; negative contribution `C=E[−Π 1_{Π<0}]=p×m` | Checked, overlap union counted once; matched-gradient association is not causality |
| Strong R90 | P90 of **positive** planar R; `A90` measured fraction of valid area; negative-transfer fraction and area-normalized enrichment `F90/A90` | Checked with valid-domain convention |
| Serial dependence | 10-frame/3000-draw moving-block-bootstrap percentile interval on frame-level medians | Checked; other paper estimands (pooled probabilities, paired null contrasts) require their own statistics, do not use this median-only helper unchanged |

### Important distinctions

1. Methods are scientifically interpretable, but dimensional limitations remain: 2D–2C observations do not yield a full 3D Liutex vector or true 3D energy-cascade flux.
2. There is no universal claim of unique or causally privileged vortex influence; gradient/strain matching indicates a shared kinematic dependence.
3. Revised Manuscript II E and SI Sec. 1.1 explicitly state U/V velocity pre-smoothing at σ=0.8 grid cells prior to separate ℓ-dependent coarse graining. This is the **default reference method**; disable it explicitly with `preSmoothSigmaCells=0` for sensitivity comparison. Exact figure-output numerical regression has not been performed.
4. Synthetic unit tests verify some mathematical identities. They cannot prove numerical measurement accuracy, full MATLAB runtime compatibility without execution, or the absence of all implementation defects.
5. The published Fig04/Fig05 producer uses a geometric centroid (MATLAB `regionprops('Centroid')`), and R-smoothing 4σ kernel support. Both are implemented in this revised reference package. Raw velocity fields and final output files remain intentionally excluded.
6. The release deliberately does NOT supply a certified full-paper figure-generation pipeline; it supplies defensible, inspectable method kernels consistent with the stated paper definitions.

### Scientific readiness verdict

**Public method-level reference source, subject to a local `run_pof('test')` pass and scientific numerical regression before claiming exact paper reproduction.** This is not an unconditional guarantee that the entire historical paper computation has zero errors. The two identified rc2 mismatches are corrected in rc3. Public open-source release separately requires permission and a license.
