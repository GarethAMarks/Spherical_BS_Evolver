# BosonStarStability

C++ code for constructing and evolving spherical boson star solutions. The core capabilities are:

- **Family construction**: sweep over central field amplitudes to build sequences of boson star models, writing ADM mass, radius, Noether charge, binding energy, frequency, and compactness to `BosonStars.dat`.
- **Radial oscillation frequencies**: linearised perturbation solver (`LinearPerturbation`) computes radial mode frequencies around a given background star.
- **Dynamical evolution in spherical symmetry**: full BSSN (or CCZ4) evolution in spherical symmetry using a modified cartoon method, with Kreiss-Oliger dissipation, fixed mesh refinement, and optional shift evolution.
- **Critical phenomena**: bisection search over initial data parameters to locate the black-hole formation threshold.

---

## Compile and run

```bash
bash make_BSEvolver_gcc.sh
./BosonStarStability.o
```

The executable must be run in the same directory as `BSParams.par`.

### Parallel runs with MPI

The time evolution can be split across MPI ranks. Build with an MPI compiler wrapper (`mpicxx`, from Open MPI or MPICH):

```bash
bash make_BSEvolver_mpi.sh
mpirun -np 4 ./BosonStarStability_mpi.o
```

This compiles the same sources with `-DUSE_MPI` into `BosonStarStability_mpi.o`, keeping its object files in `obj/MPI` so the serial and MPI builds can coexist. Without `-DUSE_MPI` every MPI call compiles to a no-op, so the serial build needs no MPI installation.

- The grid is split into contiguous blocks of points, one per rank, exchanging 3 ghost points with each neighbour per RK4 stage. Rank 0 solves for the initial star and broadcasts it, and does all file output.
- Output is the same as a serial run up to round-off: field differences are of order 1e-12, coming from the order of sums across ranks.
- Restrictions for more than one rank: the grid must be uniform (`refinement_points` disabled, i.e. first entry `<= 0`) and `adaptive_timestep = 0`, and each rank needs at least 4 points. The code stops with an error otherwise.
- Boson and Proca star evolutions both run in parallel. Only the evolution is parallelised. Use the serial build for model cycling (`cycle_only = 1`) and the linear perturbation solver (`pert_only = 1`).
- At the end of each evolution rank 0 prints the wall time of the time-stepping loop, e.g. `[timing] evolve_loop: 35.4 s over 500 timesteps (70748 us/step) on 2 ranks`, which is useful for checking scaling.

### Proca stars

Setting `proca = 1` uses the `ProcaStar` class, which reads the same model parameters as the boson star. Differences in meaning:

- `A_central` is `f(0)/alpha(0) = sqrt(-A_a Abar^a)` at the centre, for the vector potential `A = e^{-i omega t}(f dt + i g dr)`.
- `eigen` counts nodes of `g` (equivalently `E = f' - omega g`); the ground state has nodeless `g` and a single node in `f`.
- The potential family (`mu`, `lambda`, `solitonic`, `sigma`) is the same, as a function of `A_a Abar^a` in place of `|phi|^2`. With self-interactions, solutions can fail where the algebraic equation for `g` loses its root; this is reported as a warning.

Outputs are `PSdata.dat` (`r f g X phi E m`), `isotropic_proca.dat` if `isotropic = 1`, and `ProcaStars.dat` in `cycle_only` mode, with the same columns as `BosonStars.dat` plus `f0`, the value of `f(0)` with `alpha(infinity) = 1`.

A Proca star can be perturbed before evolution with `perturb = 1`, using the same `perturb_amp`, `perturb_spread` and `perturb_center` parameters as the boson star. The vector potential components `f` and `g` are multiplied by `1 + perturb_amp * exp(-(r - perturb_center)^2 / perturb_spread^2)`, then the metric and `E` are re-solved from the Hamiltonian constraint and Gauss law with `f` and `g` held fixed, so both constraints still hold. `perturb_amp > 0` makes the star denser. The perturbed profile is written to `PSdata_perturbed.dat`.

With `cycle_only = 0`, the star is evolved with the same BSSN/CCZ4 machinery, using the `ComplexVectorField` matter model (free Proca field only; `lambda != 0` or `solitonic = 1` is rejected for evolution). Evolved Proca variables, real and imaginary parts of each, are the radial components `X_z` and `E^z` (odd in `z`), `X_phi = -n^a A_a` and the Gauss-constraint damping variable `Z` (even in `z`), following Zilhão, Witek & Cardoso (2015). The modified-cartoon reduction and the 3+1 equations are checked by the scripts in `derivations/` (sympy). Not supported for Proca stars: Gaussian starts, critical studies, thin-shell data, the real scalar field, the initial Hamiltonian solver and the linear perturbation solver.

| Parameter | Default | Description |
|---|---|---|
| `proca_damping` | `1.0` | `kappa` in `d_t Z = alpha (D_i E^i + mu^2 X_phi - kappa Z) + ...` |

The Proca system is first order in space, and with centred differences it develops a grid-scale instability at the origin unless Kreiss-Oliger dissipation is on: use `damping_factor` of about 1 (0.2 is not enough). The dissipation stencil treats `c_chris_Z`, `beta`, `X_z` and `E^z` as odd across `z = 0`.

Differences in the evolution output for Proca runs:

- `SliceData.dat` and checkpoints have 8 extra trailing columns: `X_re X_im E_re E_im Xphi_re Xphi_im Z_re Z_im`.
- `Diagnostics.dat` has an extra trailing column with the magnitude of the Gauss constraint `|D_i E^i + mu^2 X_phi|`.
- In `constraint_norms.dat`, `A_central`, `phi_re_central` and the complex field energy refer to `|X_phi|`, `Re X_phi` and the Proca field. For a stationary star `|X_phi(0)| = A_central`. A trailing `Gauss_L2` column gives the Gauss constraint norm relative to the norm of `mu^2 X_phi`.

### Dirac stars

Setting `dirac = 1` (which takes priority over `proca`) uses the `DiracStar` class to construct spherical Dirac stars: two gravitating Dirac fields of mass `mu` in the `j = 1/2`, `m_j = +-1/2` s-wave states, which together give a spherically symmetric stress tensor (as in Finster, Smoller & Yau 1999 and Herdeiro, Pombo & Radu 2017). Only `D = 4` is supported. `lambda` adds a quartic (Soler-type) self-interaction and `solitonic = 1` the degenerate-vacuum potential with vacuum set by `sigma` (both described below). `cycle_only = 1` cycles a family; `cycle_only = 0` solves a single star and evolves it with the `ComplexSpinorField` matter model.

In the code's metric `ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2`, each spinor is `e^{-i omega t} alpha^{-1/2} r^{-1} (G(r) Omega_{-1,m}, i F(r) Omega_{+1,m})` and the Dirac equation reduces to

    G' = X [  G/r + (mu + omega/alpha) F ]
    F' = X [ -F/r + (mu - omega/alpha) G ]

with number density `n = (F^2 + G^2) / (2 pi alpha r^2)` and energy density `rho = omega n / alpha`; the radial pressure and the full derivation are in `derivations/check_dirac_star.py`, which also checks that the Dirac equation separates, that the stress tensor is static and spherically symmetric, and that it is conserved. Differences in parameter meaning:

- `A_central` is the central spinor amplitude `g0 = G'(0) / sqrt(alpha(0))`, which is unchanged when the lapse is rescaled to `alpha(infinity) = 1`. The maximum-mass ground state has `A_central` of about 0.23.
- `eigen` counts nodes of `G`; the ground state has nodeless `F` and `G`.
- `noether_charge` is the particle number of both spinors together.
- `lambda` is the coupling of the quartic self-interaction, described in the next subsection (`lambda = 0`: free field).
- `perturb = 1` scales both spinor amplitudes by `1 + perturb_amp * exp(-(r - perturb_center)^2 / perturb_spread^2)` and re-solves the metric (Hamiltonian constraint for `X`, polar slicing for the lapse) with the amplitudes held fixed, so the data satisfy the constraints; at `t = 0` the amplitudes are real, so the momentum density vanishes and `K_ij = 0`. `perturb_amp > 0` makes the star denser. The perturbed profile is written to `DSdata_perturbed.dat`.

Outputs are `DSdata.dat` (`r F G X phi m f g`, with `f`, `g` the physical spinor functions `F, G / (sqrt(alpha) r)`), `isotropic_dirac.dat` if `isotropic = 1`, and `DiracStars.dat` in `cycle_only` mode, with the same columns as `BosonStars.dat`. The ground-state family reaches a maximum mass `M = 0.7089 / mu` (G = 1) at `omega = 0.825 mu`, with minimum frequency `0.733 mu`; Herdeiro, Pombo & Radu (2017) give `M_max = 0.709` and minimum frequency `0.733`.

#### Quartic self-interaction

With `lambda != 0` the Lagrangian of the two spinors gains a term `-U(S)`, `U(S) = (lambda/2) S^2`, with `S = sum_m psibar_m psi_m` the scalar density (Soler model). Each spinor then obeys the Dirac equation with the effective mass `mu_eff = mu + lambda S`, so `lambda > 0` is repulsive; the stress tensor gains `g_ab (lambda/2) S^2`. In the solver variables `S = (G^2 - F^2) / (2 pi alpha r^2)`, `mu -> mu_eff` in the radial equations, and

    rho = (omega/alpha) n - (lambda/2) S^2,   p_r = [ (omega/alpha)(F^2 + G^2) - mu (G^2 - F^2) + 2 F G / r ] / (2 pi alpha r^2) - (lambda/2) S^2.

In the evolution, `S = (|u|^2 - |v|^2) / (2 pi)`, `mu -> mu_eff` in the equations for `u` and `v`, `rho = mu S + (lambda/2) S^2 + (kinetic stresses)`, and both orthonormal stresses gain `(lambda/2) S^2`. `derivations/check_dirac_star.py` checks the stationary equations and the conservation of the stress tensor with the interaction, and `derivations/check_dirac_selfint_3p1.py` checks the 3+1 sources against the evolution equations through energy and momentum conservation (both fix the coefficient of the `g_ab S^2` term).

In terms of Dzhunushaliev & Folomeev (2019), whose `lambda > 0` is attractive, `lambda_bar_DF = -lambda / (4 pi)` (`G = mu = 1`). At strong coupling (`lambda mu^2 >> 1/G`) the star approaches a fluid with `rho = mu n + lambda n^2 / 2`, `p = lambda n^2 / 2` (locally `mu_eff = omega / alpha`, `F -> 0`). Its TOV sequence has `M_max = 0.1158 sqrt(lambda)` (`G = 1`, i.e. `0.41 sqrt(|lambda_bar_DF|) M_Pl^2 / mu`) and `M / R = 0.21` at the maximum mass, stiffer than the scalar quartic (CSW) boson stars, which reach `M / R = 0.16`. With the solver (`G = mu = 1`):

| `lambda` | `M_max` | `M_max / sqrt(lambda)` | `A_central` at `M_max` | `omega` at `M_max` | `M / r_99` at `M_max` |
|---|---|---|---|---|---|
| 0 | 0.7089 | | 0.23 | 0.824 | |
| 100 | 1.309 | 0.131 | 0.17 | 0.777 | 0.159 |
| 1000 | 3.708 | 0.1172 | 0.064 | 0.755 | 0.210 |

At strong coupling the shooting is exponentially sensitive to the frequency inside the star, where small deviations of the field run away (`mu_eff` grows with the field). Diverging trial solutions are therefore stopped once `|f|` or `|g|` exceeds `3 A_central` (only for `lambda != 0`, so the free-field solver is unchanged). Even so, deep on the unstable branch the frequency cannot be resolved in double precision and the solver fails (for `lambda = 1000`, beyond `A_central` of about 0.11, about twice the maximum-mass amplitude).

Evolution tests at `lambda = 1000` (`R = 150`, `dr = 0.05`, `BS_resolution_factor = 4`, `isotropic = 1`, `damping_factor = 1`), on either side of the maximum mass at `A_central = 0.064`:

- Stable branch, `A_central = 0.04` (`M = 3.216`): unperturbed, `|u(0)|` stays constant to `1e-5` up to `t = 1000`. With `perturb_amp = +-0.02` the star oscillates about equilibrium with a non-growing amplitude (about 1% of `A_central`).
- Unstable branch, `A_central = 0.09` (`M = 3.420`):
  - `perturb_amp = +0.01` collapses, with an apparent horizon at `t = 121.5` (`121.2` at `dr = 0.025`).
  - Unperturbed, truncation error triggers collapse at `t = 275.5`.
  - `perturb_amp = -0.01` migrates: the star expands to `|u(0)| = 0.025` and then oscillates about the stable-branch star with the same particle number (`A_central = 0.042`). The particle number is conserved to `1e-5` until ejected matter reaches the outer boundary. The trajectory agrees with `dr = 0.025` to `1e-5`.

#### Degenerate-vacuum potential

`solitonic = 1` gives the Dirac analogue of the solitonic boson-star potential. The mass term and self-interaction combine into

    W(S) = mu S + U(S) = mu S (1 - S/S0)^2,     S0 = sigma^2 / (4 pi),

which vanishes at `S = 0` and at the second (degenerate) vacuum `S = S0`. At the centre of a star `S = A_central^2 / (2 pi)`, so the vacuum is at `A_central = sigma / sqrt(2)`, as for the solitonic boson stars. The fermions are massless in that vacuum: `mu_eff = mu (1 - S/S0)(1 - 3 S/S0)`. In the notation of Herdeiro & Radu (2020), this is their quartic-sextic potential with `beta = 1/4`. It can be combined with `lambda`. The potential is implemented in `DiracPotential.h` and used by both the solver and the evolution. The symbolic checks cover any `U = (lambda/2) S^2 + kappa S^3`, which includes this case.

The shooting solver works for `sigma` of about 0.3 and above. For smaller `sigma` the stars develop thin walls and it fails. `relaxation/` contains a relaxation solver for these models:

- `dirac_relax.py` solves the same equations, written for the physical amplitudes `f`, `g`, with `m`, `phi` and the frequency promoted to a grid variable (`omega' = 0`).
  - Boundary conditions: regularity at the centre, `alpha -> 1`, and exponential decay.
  - Method: scipy's `solve_bvp` (4th-order collocation). The mesh is redistributed at every step so that thin walls stay resolved.
  - Families are followed by pseudo-arclength continuation in `(A_central, omega, ln M)`.
- `families.py <sigma> <dir> <binary>` computes a whole family. It is seeded by one shooting solution at `A_central = sigma / (2 sqrt 2)`. It writes `DiracStars_sigma<sigma>.dat`, ordered along the family:
  - first the columns of `BosonStars.dat`: `A_central M r_99 noether_charge binding_energy omega omega_unrescaled phi0 C4 Cmax`;
  - then `unstable_modes` (from the turning-point analysis), `n_light_rings` and two energy-condition measures, `wec = min rho / max rho` and `dec = max(|p_r|, |p_t|) - rho` over `max rho`.
  - Without shooting, `omega_unrescaled` has no unique meaning, since it depends on the initial guess for the central lapse. It is given as `omega / alpha(0)`.
- `export_profile.py <dir> <sigma> list | <index> <file>` lists the stored solutions of a family (every 10th point, with their stability and light rings) or exports one for evolution. The exported file has `#` header lines with `omega`, `sigma`, `lambda`, `mu` and `solitonic`, then the columns `r f g m phi` on a uniform grid. It is evaluated from the collocation interpolant after re-solving at tolerance `1e-8`. Setting `dirac_profile = 1` and `dirac_profile_file = <file>` in `BSParams.par` makes `DiracStar` read this profile instead of shooting:
  - It is interpolated onto the evolution grid with 4-point Lagrange polynomials, with a Schwarzschild exterior beyond its last radius.
  - The potential parameters must match `BSParams.par`.
  - Everything downstream works as for a shot star: `BS_resolution_factor`, isotropic coordinates, `perturb` and evolution.
- `stability.py` applies the turning-point method. Stability can only change at the joint extrema of `M` and `N`, which are cusps in the `(N, M)` plane. At each cusp, the branch with the larger `M - omega_c N` at a given `N` has one more unstable mode. The count starts from the stable Newtonian end.
- `plot_families.py` plots the compactness and mass along the families and marks:
  - the ends of the stable branches;
  - the onset of light rings;
  - where the energy conditions first fail.

For `sigma = 0.3` the relaxation solver agrees with shooting to a few parts in `1e6` in `M`, `N` and `omega`. Halving the tolerance, moving the inner boundary or enlarging the domain changes `M` by about `1e-8`. Along the families `dM = omega dN` holds to about 1% (the residual is from finite differences between neighbouring models).

Each family has the following structure:
- a tiny Newtonian stable branch, ending at a maximum mass of about `0.5 sigma` (for `sigma <= 0.07` it lies below the computed amplitudes);
- an unstable branch, down to a minimum of the mass at `A_central = 0.43 sigma/sqrt 2`;
- a long stable branch from that minimum to the main maximum mass, slightly beyond the vacuum amplitude;
- an unstable spiral after the main maximum.

Results along the main stable branch (`G = mu = 1`, `C = M / r_99`, `C_max = max m(r)/r`):

| `sigma` | `M_max` | `M_max sigma^2` | `omega` | `C` | `C_max` | `A_central / (sigma/sqrt 2)` | `C` at light-ring onset |
|---|---|---|---|---|---|---|---|
| 0.3 | 2.810 | 0.253 | 0.1101 | 0.416 | 0.417 | 1.277 | 0.322 |
| 0.2 | 5.812 | 0.233 | 0.0634 | 0.434 | 0.430 | 1.213 | 0.333 |
| 0.15 | 9.655 | 0.217 | 0.0432 | 0.442 | 0.437 | 1.188 | 0.337 |
| 0.1 | 19.63 | 0.196 | 0.0252 | 0.448 | 0.444 | 1.147 | 0.340 |
| 0.07 | 36.54 | 0.179 | 0.0156 | 0.451 | 0.448 | 1.116 | 0.341 |
| 0.05 | 65.64 | 0.164 | 0.0099 | 0.452 | 0.451 | 1.099 | 0.336 |
| 0.035 | 122.2 | 0.150 | 0.0060 | 0.454 | 0.453 | 1.082 | 0.341 |
| 0.025 | 219.8 | 0.137 | 0.0038 | 0.454 | 0.454 | 1.069 | 0.340 |
| 0.018 | > 388.2 | > 0.126 | < 0.0024 | > 0.453 | > 0.453 | > 1.033 | 0.338 |
| 0.0125 | > 696.9 | > 0.109 | < 0.0015 | > 0.445 | > 0.445 | > 1.009 | 0.337 |

The maximum mass ends the stable branch. For `sigma = 0.018` and `0.0125` the continuation stalls on the stable branch just before the maximum, so those entries are bounds. These runs reached the following (on two cores):
- `sigma = 0.025`: about 35 minutes for the whole family.
- `sigma = 0.0125`: `M = 697` and `C = 0.445`, with the walls about `1e-3` of the stellar radius.

The table and the plot come from

    python3 relaxation/families.py 0.1 <dir> ./BosonStarStability.o
    python3 relaxation/plot_families.py <dir> families.png

Light rings appear on the stable branch at `C` of about 0.34, so a large part of every stable branch has light rings. These stars are expected to suffer the nonlinear light-ring instability, which the turning-point method does not capture.

These compactnesses exceed the Buchdahl bound 4/9, which assumes isotropic, positive pressure no larger than the density. The classical Dirac energy density `rho = (omega/alpha) n - P(S)` is not positive definite, and the cores of these stars violate the energy conditions:
- **Dominant energy condition:** `|p| <= rho` first fails at `omega` of about 0.35 for every `sigma`, at `C` from 0.16 (`sigma = 0.3`) down to 0.0005 (`sigma = 0.0125`). Beyond that, the pressure in the core exceeds the energy density by large factors.
- **Weak energy condition:** the energy density becomes negative in the core. The minimum, relative to the maximum of `rho`, is about -2% for `sigma = 0.3` and about -3e-4 for `sigma = 0.0125`. Beyond 0.1% it first appears at `C` of about 0.13–0.26 for `0.035 <= sigma <= 0.3`.
- **Thin-wall limit:** in the locally homogeneous approximation the interior sits just above the vacuum, with `mu_eff = omega/alpha`. There `rho = -2 mu S0 (S/S0 - 1)^2 < 0` and `p = omega S0 / alpha > 0`, so the mass sits in the wall.

These are classical-field solutions. For fermions, a quantised treatment would not allow negative energy densities.

In the evolution the tetrad is fixed by the time gauge (`e_0 = n`) and the symmetric spatial triad, which in spherical symmetry is aligned with the diagonal `h_zz`, `h_ww`, so it is a function of the BSSN variables and no extra variables are evolved. The evolved variables are the complex orthonormal-frame amplitudes `u` (even in `z`) and `v` (odd), with `u = e^{-i omega t} g`, `v = e^{-i omega t} f` for a stationary star:

    d_t u = beta u' + (alpha K/2 - i alpha mu) u + i alpha [ sqrt(chi/h_zz) (v' + L v) + (sqrt(chi/h_zz) + sqrt(chi/h_ww)) v / z ]
    d_t v = beta v' + (alpha K/2 + i alpha mu) v - i alpha [ sqrt(chi/h_zz) (u' + L u) + (sqrt(chi/h_zz) - sqrt(chi/h_ww)) u / z ]

with `L = (h_ww'/h_ww - chi'/chi + alpha'/alpha) / 2` (and `mu -> mu_eff` with the self-interaction). These and the matter sources are derived and checked in `derivations/check_dirac_3p1.py`. Like the Proca system the equations are first order in space, and the evolution is unstable at the origin without Kreiss-Oliger dissipation: use `damping_factor` of about 1. MPI runs are supported (uniform grid, as for the other matter models).

Differences in the evolution output for Dirac runs:

- `SliceData.dat` and checkpoints have 4 extra trailing columns: `u_re u_im v_re v_im`.
- In `constraint_norms.dat`, `A_central`, `phi_re_central` and the complex field energy refer to `|u(0)|`, `Re u(0)` and the Dirac field, and the Noether charge is the particle number of both spinors. For a stationary star `|u(0)| = A_central`.

---

## Output files

| File | Written every | Contents |
|---|---|---|
| `BSdata.dat` | once, at startup | Radial profiles from the ODE solver: `r A X phi eta m` |
| `isotropic.dat` | once, if `isotropic = 1` | Isotropic-coordinate field values |
| `SliceData.dat` | `write_interval` steps | Radial profiles of all 13 BSSN + matter variables (see `bssn_var` enum in `EvolutionVariables.h`) |
| `Diagnostics.dat` | `write_interval` steps | Radial profiles of Hamiltonian constraint, momentum constraint, conformal-metric determinant, auxiliary variable, field amplitude, geometric and matter Ricci scalars |
| `constraint_norms.dat` | `write_CN_interval` steps | Time-series: L² constraint norms, central chi, central field amplitude, AH radius, ADM mass, Noether charge, central dK/dt, central frequency (from the phase change of the central field over one step), central lapse, field energies, central 4-Ricci scalar, central Re(field), AH mass |
| `BosonStars.dat` | once per model, in `cycle_models` | Family data: `A_central M r_99 noether_charge binding_energy omega omega_unrescaled phi0 C4 Cmax` |
| `checkpoint<T>.dat` | when `t` crosses a multiple of `checkpoint_time` | Full slice snapshot for restarting |
| `Diagnostics_<T>.dat`, `SliceData_<T>.dat` | when `t` crosses a multiple of `dump_time` | Snapshot diagnostics/slice at time `T` |

---

## BSParams.par reference

### Boson star model

| Parameter | Default | Description |
|---|---|---|
| `proca` | `0` | Construct (and evolve) Proca stars instead of boson stars; see below |
| `D` | `4` | Spacetime dimension |
| `mu` | `1.0` | Scalar field mass |
| `lambda` | `0` | Quartic self-coupling coefficient (Dirac stars: `U = (lambda/2) (psibar psi)^2`, see above) |
| `solitonic` | `0` | Use solitonic potential `V = mu^2 (1 - 2(A/sigma)^2)^2`; if false, use mini-BS potential. Dirac stars: the degenerate-vacuum potential `mu S (1 - S/S0)^2`, `S0 = sigma^2 / (4 pi)` (see Dirac stars) |
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
| `adaptive_timestep` | `0` | With refinement, subcycle the refinement levels in time (Berger-Oliger): each level steps at `courant_factor` times its own spacing, so coarse levels take larger steps. See below |

**Subcycling (`adaptive_timestep = 1`).** The coarsest level takes one step of `courant_factor * dr * 2^(L-1)` (for `L` levels), then each finer level takes two steps of half the size, recursively. Each level also evolves a 14-point buffer (at its own spacing) into its neighbours. Before each substep, the buffer reaching into the next coarser level is filled using RK4 dense output in time and 6-point Lagrange interpolation in space. Notes:
- Requires refinement and a serial run. Each refinement level must be wide enough to hold the buffer of the next finer level (the code stops with an error otherwise).
- `write_CN_interval`, `write_interval` and `AH_search_steps` count coarse steps.
- For refined runs with mostly coarse grid, it is typically 1.6-2.6x faster at the same or better accuracy than `adaptive_timestep = 0`.
- With `adaptive_timestep = 0` (the default) results are unchanged.

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
| `damping_factor` | `0.0` | Kreiss-Oliger dissipation strength (about 1 is needed for Proca stars) |
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
