# BosonStarStability

C++ code for constructing and evolving spherical boson star solutions. The core capabilities are:

- **Family construction**: sweep over central field amplitudes to build sequences of boson star models, writing ADM mass, radius, Noether charge, binding energy, frequency, and compactness to `BosonStars.dat`.
- **Radial oscillation frequencies**: linearised perturbation solver (`LinearPerturbation`) computes radial mode frequencies around a given background star.
- **Dynamical evolution in spherical symmetry**: full BSSN (or CCZ4) evolution in spherical symmetry using a modified cartoon method, with Kreiss-Oliger dissipation, adaptive mesh refinement via manual refinement points, and optional shift evolution.
- **Critical phenomena**: bisection search over initial data parameters to locate the black-hole formation threshold.

---

## Compile and run

```bash
bash make_BSEvolver_gcc.sh
./BosonStarStability
```

The executable must be run in the same directory as `BSParams.par`.

---

## Output files

| File | Written every | Contents |
|---|---|---|
| `BSdata.dat` | once, at startup | Radial profiles from the ODE solver: `r A X phi eta m` |
| `isotropic.dat` | once, if `isotropic = 1` | Isotropic-coordinate field values |
| `SliceData.dat` | `write_interval` steps | Radial profiles of all 13 BSSN + matter variables (see `bssn_var` enum in `EvolutionVariables.h`) |
| `Diagnostics.dat` | `write_interval` steps | Radial profiles of Hamiltonian constraint, momentum constraint, conformal-metric determinant, auxiliary variable, field amplitude, geometric and matter Ricci scalars |
| `constraint_norms.dat` | `write_CN_interval` steps | Time-series: L² constraint norms, central chi, central field amplitude, AH radius, ADM mass, Noether charge, field energies, central 4-Ricci scalar, AH mass |
| `BosonStars.dat` | once per model, in `cycle_models` | Family data: `A_central M r_99 noether_charge binding_energy omega omega_unrescaled phi0 C4 Cmax` |
| `checkpoint<T>.dat` | when `t` crosses a multiple of `checkpoint_time` | Full slice snapshot for restarting |
| `Diagnostics_<T>.dat`, `SliceData_<T>.dat` | when `t` crosses a multiple of `dump_time` | Snapshot diagnostics/slice at time `T` |

---

## BSParams.par reference

### Boson star model

| Parameter | Default | Description |
|---|---|---|
| `D` | `4` | Spacetime dimension |
| `mu` | `1.0` | Scalar field mass |
| `lambda` | `0` | Quartic self-coupling coefficient |
| `solitonic` | `0` | Use solitonic potential `V = mu^2 (1 - 2(A/sigma)^2)^2`; if false, use mini-BS potential |
| `sigma` | `0.2` | Solitonic potential parameter (only used when `solitonic = 1`) |
| `G` | `1.0` | Newton's constant |
| `eigen` | `0` | Eigenstate index: 0 = ground state, 1 = first excited state, etc. |
| `alpha_central` | `0.5` | Initial central lapse for the ODE solver (rescaled after shooting) |
| `A_central` | `0.05` | Central field amplitude for the boson star |
| `frequency_guess` | `1.8` | Initial frequency guess for the shooting method |
| `freq_epsilon` | `1e-19` | Convergence tolerance for the frequency |
| `isotropic` | `0` | If true, convert initial data to isotropic coordinates before evolution |
| `r_match_fac` | `0.85` | Fraction of the grid outer boundary at which the BS exterior solution is matched |

### Dynamical perturbation

| Parameter | Default | Description |
|---|---|---|
| `perturb` | `0` | Add a Gaussian perturbation to the field before evolution |
| `perturb_amp` | `0.05` | Amplitude of the Gaussian perturbation |
| `perturb_spread` | `6.0` | Width (sigma) of the Gaussian perturbation |
| `perturb_center` | `0` | Radial centre of the Gaussian perturbation |
| `mirror_gaussian` | `0` | Mirror the Gaussian about `r = 0` for a purely ingoing pulse |
| `gaussian_start` | `0` | Replace the BS initial data entirely with a pure Gaussian scalar-field pulse |
| `enforced_freq` | `0.9` | Oscillation frequency used when `gaussian_start = 1` |

### Grid and timestepping

| Parameter | Default | Description |
|---|---|---|
| `R` | `200` | Outer boundary radius |
| `n_gridpoints` | `2000` | Number of radial grid points |
| `courant_factor` | `0.4` | CFL factor: `dt = courant_factor * dr` |
| `stop_time` | `1000` | Coordinate time at which to stop the evolution |
| `cell_centered` | `0` | If true, use a cell-centred grid (first point at `dr/2`); otherwise vertex-centred |
| `BS_resolution_factor` | `1` | Multiply the BS ODE solver resolution by this factor (must be a power of 2) |
| `refinement_points` | `(none)` | Space-separated list of grid indices after which the effective resolution is halved; set first entry negative to disable refinement |
| `adaptive_timestep` | `0` | Reduce the timestep at higher refinement levels |

### Massless real scalar field

| Parameter | Default | Description |
|---|---|---|
| `add_real_field` | `0` | Inject an in-falling Gaussian real scalar field pulse before evolution begins |
| `real_amp` | `0.01` | Amplitude of the real scalar field pulse |
| `real_sigma` | `5.0` | Width of the real scalar field pulse |
| `real_center` | `50.0` | Radial centre of the real scalar field pulse |

### Critical phenomena

| Parameter | Default | Description |
|---|---|---|
| `do_ah_search` | `1` | Search for an apparent horizon at every timestep |
| `critical_study` | `0` | Enable bisection search for the black-hole formation threshold |
| `hi_guess` | `0.05` | Initial supercritical guess for the bisection parameter |
| `lo_guess` | `0.001` | Initial subcritical guess for the bisection parameter |
| `critical_eps` | `1e-12` | Bisection terminates when `hi - lo < critical_eps` |
| `lapse_thresh` | `(large)` | (Unused by default) Central lapse value below which the run is declared supercritical |
| `A_thresh` | `(large)` | Central field amplitude above which the run is declared supercritical (alternative to AH search) |
| `sub_min_time` | `75` | Minimum coordinate time before a run can be declared subcritical |
| `subcritical_time` | `1000` | If no horizon has formed by this time, the run is declared subcritical |
| `fill_subcritical` | `0` | After bisection completes, run additional subcritical evolutions to fill out the phase diagram |

### Evolution

| Parameter | Default | Description |
|---|---|---|
| `run_spacetime_solver` | `0` | Re-solve the spacetime constraints at every step (testing only) |
| `min_chi` | `0.0001` | Floor on the conformal factor chi to prevent division by zero near collapse |
| `min_z` | `0.00001` | Below this radius, regularised (`z → 0`) forms of the equations are used |
| `max_stored_slices` | `5` | Number of time slices held in memory simultaneously (minimum 3 for RK4) |
| `stop_on_migrate` | `0` | Exit with a message if the central field amplitude changes by more than 10% |
| `sigma_BSSN` | `0.333333` | Coefficient of the constraint-damping term in the conformal connection evolution; in CCZ4 mode, also sets `c3 / 3` |
| `damping_factor` | `0.0` | Kreiss-Oliger dissipation strength |
| `only_BS_violation` | `1` | Restrict constraint damping to the interior of the star (`r < r_99`) |
| `spatially_varying_BC` | `1` | Use spatially varying (Sommerfeld-like) outer boundary conditions |
| `shock_gauge` | `0` | Use shock-avoiding Bona-Masso lapse with `f(alpha) = 1 + shock_fac / alpha^2` instead of 1+log |
| `shock_fac` | `2.0` | Parameter in the shock-avoiding gauge condition (only used when `shock_gauge = 1`) |
| `one_log_fac` | `2.0` | Coefficient in the 1+log slicing condition `d_t alpha = -one_log_fac * alpha * K + ...` |
| `eta` | `1.0` | Gamma-driver damping coefficient in the shift evolution equation |
| `gamma_fac` | `0.75` | Coefficient of the conformal connection in the Gamma-driver shift condition |
| `evolve_shift` | `0` | If false, the shift is fixed to zero throughout |
| `use_CCZ4` | `0` | Use CCZ4 formulation instead of BSSN |
| `c1` | `0.02` | CCZ4 constraint-damping coefficient `kappa_1` |
| `c2` | `0.0` | CCZ4 constraint-damping coefficient `kappa_2` |
| `theta_off` | `0` | Disable evolution of the CCZ4 theta variable (keep off for standard use) |

### Output and checkpointing

| Parameter | Default | Description |
|---|---|---|
| `start_time` | `0` | Coordinate time label assigned to the initial slice (useful when restarting) |
| `write_interval` | `1` | Write `SliceData.dat` and `Diagnostics.dat` every this many timesteps |
| `write_CN_interval` | `5` | Append a row to `constraint_norms.dat` every this many timesteps |
| `checkpoint_time` | `2000` | Write a full checkpoint slice whenever `t` passes a new multiple of this value |
| `dump_time` | `0` | Write timestamped `SliceData_<T>.dat` and `Diagnostics_<T>.dat` snapshots whenever `t` passes a new multiple of this value; disabled if zero |
| `run_quietly` | `0` | Suppress per-step console output |
| `store_A0` | `1` | Store the initial central field amplitude for later use in migration detection |
| `read_thinshell` | `0` | Read initial data from the output of an external thin-shell solver instead of the built-in BS ODE solver |
| `cutoff_frac` | `0.75` | Fraction of the grid used when reading thin-shell initial data |
| `make_tangherlini` | `0` | Initialise a Tangherlini (higher-dimensional Schwarzschild) black hole instead of a boson star |
| `uniform_data` | `1` | Expect thin-shell data pre-interpolated to a uniform grid via `smooth_and_interp.py` |

### Linear perturbation solver

| Parameter | Default | Description |
|---|---|---|
| `pert_only` | `0` | Run the linear perturbation solver only (skip dynamical evolution) |
| `chi_sq0` | `0.0` | Initial value of chi² for the perturbation ODE |
| `gamma0` | `0.0` | Initial frequency-squared guess for the perturbation eigenvalue |
| `chi_range` | `0.15` | Range over which to scan chi² for the perturbation modes |
| `chi_epsilon` | `1e-8` | Convergence tolerance for the perturbation shooting |
| `cutoff_radius` | `30` | Radius beyond which the perturbation solution is matched to the exterior |
| `mode_number` | `2` | Number of nodes of the perturbation eigenfunction to target |
| `noether_epsilon` | `1e-5` | Tolerance for Noether charge normalisation of perturbation modes |
| `dynamic_cutoff` | `0` | Dynamically adjust the cutoff radius during the perturbation solve |

### Model cycling

| Parameter | Default | Description |
|---|---|---|
| `cycle_only` | `1` | Only construct the family of models without running a dynamical evolution or perturbation calculation |
| `A0` | `0.001` | Starting central amplitude for the model family sweep |
| `dA` | `0.0005` | Step size in central amplitude between successive models |
| `n_stars` | `1800` | Number of models to construct in the sweep |
