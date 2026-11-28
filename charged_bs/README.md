# Static charged boson stars → SBSE initial data

Static charged solitonic boson stars (Einstein-Maxwell-Klein-Gordon, thin-shell regime) and the tools that turn them
into initial data for `charged = 1` runs. Units are G = μ = 1, Gaussian Maxwell, V = |φ|²(1 − 2|φ|²/σ0²)², and q = q̂ σ0.

| File | Role |
|---|---|
| `emkg.py` | Static equations, observables and grids. Uses D = ∇ − i q A and φ = A e^{−iωt}. |
| `emkg2.py` | Moving-mesh Newton solver, with the wall pinned at a node. |
| `nat2.py`, `cbs_lo.py` | Leading-order (σ0 → 0) seed for the first solution. |
| `arc2.py` | Pseudo-arclength continuation in Wc = Ω_c / σ0. Has an extra `snap` callback that keeps the profiles. |
| `make_charged_id.py` | `trace` / `pairs` / `write`: family tracing, same-N pairs, and `charged_BSdata.dat` for SBSE. |
| `plot_convergence.py` | Three-resolution constraint convergence plots from `constraint_norms.dat`. |

The solver files are copied from the charged thin-shell project: its bundle, with `arc2.py` from the same-N pair study.

**Sign convention.** SBSE's stars rotate as φ ∝ e^{+iωt}, the complex conjugate of `emkg.py`'s convention. So SBSE's
`ChargedScalarField` uses D = ∇ + i q A. The gauge field A_t = Φ(r) is the same in both, and q > 0 gives Q > 0. The file
written by `make_charged_id.py` holds the gauge-invariant data A, X, ln α, E^z = Q/(r² X), Φ_e = −A_t/α and
Ω = (ω + q A_t)/α, with A_t → 0 and α → 1 at infinity. SBSE sets K_φ = −½ i Ω A.

**Stability along a family.** The turning-point rule applies: one radial mode changes sign at each extremum of N along the
family. On the families traced so far, the side of the N maximum with the larger central lapse is the stable one.
