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

**Evolutions.** Models read from relaxation profiles (`dirac_profile = 1`) were evolved with `isotropic = 1` and `damping_factor = 1`. They start stationary, with the Hamiltonian constraint at about `1e-11` inside the star. The turning-point method does not capture their dynamics: much of the turning-point-stable branch is dynamically unstable.

- **Unstable modes:** the instability is oscillatory, a complex eigenvalue as known for solitary waves of nonlinear Dirac equations. It does not appear at an extremum of `M`, so the turning-point method cannot see it.
- **Resolution checks:**
  - It is seeded by truncation error. For `sigma = 0.3`, `C = 0.137`, the growth rate is 0.307, 0.321 and 0.320 at `dr = 0.05, 0.025, 0.0125`, with period about 8.
  - The growing perturbation satisfies the momentum constraint to 4th order. It is a physical solution, not a constraint-violating mode.
- **`sigma = 0.3`:** models with `C = 0.014, 0.029, 0.052` stay stationary. The deviation stays at the truncation level up to `t = 200`, and to `t = 2000` for `C = 0.052`. At `C = 0.086` (`A_central = 0.91 A_vac`) a slow instability appears, and at `C = 0.137` it e-folds in about 3.
- **`sigma = 0.07`, light-ring regime:** all three models collapse to black holes. The apparent horizon forms with essentially the whole mass, enclosing the frozen core (central lapse about `0.03`).
  - The perturbation grows just inside the wall. Its growth rate rises from about 0.05 to 0.20 by `t = 60`, as the fastest mode takes over.
  - `C = 0.403` (light rings): horizon at `t = 65`.
  - `C = 0.450` (Buchdahl-violating, turning-point stable): horizon at `t = 99.5`, or `104.3` at `dr = 0.025`.
  - `C = 0.449` (just past the maximum mass): horizon at `t = 108`.

**Radial pulsation modes.** `derivations/dirac_pulsation.py` derives the linearised radial perturbation equations for two spinors with any `U(S)`, in the polar-areal gauge `ds^2 = -alpha^2 dt^2 + X^2 dr^2 + r^2 dOmega^2`.
- It computes the Einstein tensor and checks the Hamiltonian constraint, the momentum constraint `X_t = -4 pi r alpha X^2 j` and the polar slicing condition.
- It linearises the Dirac equations in the frame rotating with the star, `u = e^{-i omega t} (g + uR + i uI)` and `v = e^{-i omega t} (f + vR + i vI)`, together with `dX` and `dphi`. All perturbations go as `e^{lambda t}`, with `X_t` eliminated by the momentum constraint.
- It checks that the background satisfies the `DiracStar` equations and that the global phase rotation is a zero mode.
- It writes the coefficients to `relaxation/dirac_pulsation_coeffs.py`.
- With `lambda = -i sigma`, real `sigma` is a stable mode: the field components oscillate at `omega +- sigma`. `Re lambda > 0` is unstable, and the instability is oscillatory when `Im lambda != 0`.

`relaxation/dirac_modes.py` (`ModeSolver`) solves the generalised eigenproblem `A Y = lambda B Y`:
- The six fields are discretised by Chebyshev collocation on `[-R, R]`, folded by parity (`u`, `dX`, `dphi` even; `v` odd).
- Boundary conditions: `u = v = 0` at `R`, `dX(0) = 0` and `dphi(R) + dX(R)/X(R) = 0`.
- The collocation coordinate mixes `r` with the tortoise coordinate `int X/alpha dr` (weight 0.8). Without this, the short local wavelengths in the redshifted cores of compact stars are unresolved and give spurious unstable eigenvalues.

`relaxation/scan_modes.py <dir> <sigma> <out.npz> [stride] [a_lo] [a_hi] [N1] [N2]` computes the spectrum along a family from `families.py`, which must be run with a small `profile_stride` (6th argument) so that the profiles are stored.
- It keeps eigenvalues that agree between two resolutions, and frequencies the coarser grid resolves.
- It records the largest real and oscillatory growth rates and the three lowest stable frequencies.

`relaxation/shoot_modes.py <binary> <out.npy> <A1,A2,...> [lambda] [solitonic] [sigma]` does the same for single stars found by shooting (free, quartic, `sigma >~ 0.3`).

Weakly unstable modes at high frequency (`Re lambda ~ 1e-3`, `|Im lambda| ~ 2`) come from the reflecting outer boundary: their growth rate falls as `R` grows. Only growth rates above 0.01 are counted below, except where an onset is extrapolated.

Validation against the 1D evolutions:

| model | eigenvalue | evolution |
|---|---|---|
| quartic `lambda = 1000`, `A = 0.04`: `sigma_1, sigma_2, sigma_3` | 0.0184, 0.0474, 0.0665 | 0.0188, 0.0474, 0.067 (FFT of `A_central(t)`) |
| quartic `A = 0.09` (past the maximum mass) | real 0.0330 | growth 0.034–0.036 |
| `sigma = 0.3`, `C = 0.086` | 0.037 + 1.678 i | envelope growth 0.0375 |
| `sigma = 0.3`, `C = 0.137` | 0.375 + 0.769 i | 0.379–0.386 + 0.761–0.763 i (fit of `|u(0)|` in the linear regime, `dr = 0.0125`) |
| `sigma = 0.07`, `C = 0.450` | real 0.199, localised just inside the wall | late-time growth 0.199, near the wall |

Results (`G = mu = 1`):
- **Free and quartic stars:**
  - The fundamental `sigma_1^2` decreases through zero at the maximum mass (`A_central = 0.23` free; about `0.063` for `lambda = 1000`) and becomes a real unstable mode, so the turning-point method is right for these families.
  - Free stars: `sigma_1` = 0.0078, 0.0258, 0.0340, 0.0269, 0.0121 at `A_central` = 0.02, 0.08, 0.14, 0.20, 0.225.
  - Free stars past the maximum: `lambda` = 0.020 at `A_central = 0.24` and 0.090 at 0.35.
  - No oscillatory instabilities.
- **Degenerate potential:** the turning-point-stable branch has two dynamical instabilities that the turning-point method misses.
  - **Oscillatory instability:** a complex quartet `lambda = +-Re +- i Im` sets in at `A_central` of about `0.89 A_vac` and `omega` of about 0.52–0.56 for every `sigma`. It leaves the high-frequency part of the spectrum (`Im lambda` about 1.8–2.1 at onset), with `Re lambda` growing linearly along the family. Its growth rate peaks near `A_central = A_vac` and falls towards the maximum mass. This happens long before light rings form. The onsets below come from spectra at every stored model, extrapolating `Re lambda` to zero:

    | `sigma` | `omega` at onset | `M` at onset | `C` at onset | `Im lambda` at onset | peak `Re lambda` | `C` at the peak |
    |---|---|---|---|---|---|---|
    | 0.3 | 0.522 | 0.341 | 0.080 | 1.82 | 0.47 | 0.17–0.19 |
    | 0.2 | 0.542 | 0.177 | 0.038 | 1.93 | 0.55 | 0.1 |
    | 0.1 | 0.557 | 0.048 | 0.010 | 2.05 | 0.61 | 0.03 |
    | 0.07 | 0.556 | 0.024 | 0.005 | 2.01 | 0.63 | 0.02 |

  - **Real (non-oscillatory) unstable modes** at larger `C`, localised at the inner edge of the wall (none for `sigma = 0.3`). The lowest-frequency complex quartet reaches the real axis (`Im lambda -> 0`) and splits into two real pairs. The larger one decreases towards the maximum mass, to 0.18 at `C = 0.447` for `sigma = 0.1` and 0.20 at `C = 0.448` for `sigma = 0.07`. For `sigma = 0.2` the two real pairs merge again into a complex quartet at `C = 0.345`.

    | `sigma` | `omega` | `M` | `C` | `lambda` when the quartet reaches the real axis |
    |---|---|---|---|---|
    | 0.2 | 0.209 | 1.336 | 0.196 | 0.43 |
    | 0.1 | 0.264 | 0.347 | 0.049 | 0.58 |
    | 0.07 | 0.269 | 0.178 | 0.025 | 0.60 |

  - The number of strongly unstable modes grows along the branch, from one at the onset to 5–10 with `|lambda| < 2.5` at `C` of about 0.1. At least some enter from frequencies above the scan window (`|lambda| = 2.5`), so their first appearance is not located here. Near the maximum mass the deeply redshifted core supports a dense band of low-frequency unstable modes, with `Re lambda` up to about 0.1.
  - Beyond `Im lambda` of about 3, a weak band of evenly spaced unstable modes (spacing about `pi / R`) also appears once the core is unstable. Its growth rate depends on the domain size, so it is not an eigenmode of the infinite problem. It may reflect convective amplification in the core.
  - The small Newtonian stable branches (`omega -> mu`) have no resolved unstable modes.
  - Only the low-compactness end of each stable branch is dynamically stable: `C < 0.080` for `sigma = 0.3`, and `C < 0.005` for `sigma = 0.07`. All the light-ring and Buchdahl-violating models are linearly unstable, with e-folding times between about 3 and 8 (in units of `1/mu`).
  - The lowest stable frequency goes to zero at both turning points. Between them, for most of the low-compactness part of each branch, it sits at the continuum edge `mu - omega`: it is then a box mode of the discretised continuum, not a bound internal mode. A bound mode below the edge separates at larger `C` (about 0.04–0.06 for `sigma = 0.3`, 0.02–0.03 for `sigma = 0.07`), reaching 0.55–0.77 before it falls to zero at the maximum mass.
  - The continuous spectrum is `|Im lambda| >= mu - omega`, with embedded thresholds at `+-(mu + omega)`: one or both sidebands `omega +- sigma` of a mode radiate. Near the onset of the oscillatory instability, the growth rate depends on the domain size by about 20% (for example 0.034 at `R = 34` and 0.037 at `R = 69` for `sigma = 0.3`, `C = 0.086`). The onset is therefore located to about `+-0.01` in `A_central / A_vac`. The evolution at `R = 160` grows at 0.0375.

**Free stars: radial modes of the first excited state.** With `eigen = 1` (one node in `G`; built by shooting and re-solved with `dirac_relax.py`), the radial spectrum has an oscillatory unstable quartet at every amplitude, the analogue of the decay of excited boson stars towards the ground state. A real unstable mode appears in addition past the excited family's maximum mass (`M = 1.410` at `A_central` about 0.2).

| `A_central` | `omega` | `M` | oscillatory `lambda` | real `lambda` |
|---|---|---|---|---|
| 0.02 | 0.9822 | 0.729 | 0.0018 + 0.0109 i | – |
| 0.05 | 0.9575 | 1.059 | 0.0048 + 0.0277 i | – |
| 0.11 | 0.9153 | 1.329 | 0.0115 + 0.0570 i | – |
| 0.17 | 0.8821 | 1.406 | 0.0189 + 0.0817 i | – |
| 0.20 | 0.8687 | 1.410 | 0.0227 + 0.0923 i | 0.0062 |
| 0.30 | 0.8373 | 1.348 | 0.0354 + 0.1189 i | 0.0323 |
| 0.40 | 0.8226 | 1.245 | 0.0467 + 0.1320 i | 0.0583 |

### Non-radial (axisymmetric) perturbations of free Dirac stars

`derivations/dirac_nonradial.py` derives the linearised equations for perturbations of total angular momentum `L` beyond spherical symmetry (free field, `D = 4`):
- **Symmetry.** The background pair (`m_j = +-1/2`) is invariant under rotations combined with a flavour rotation of the doublet, so modes carry an integer `L`. In the axisymmetric (`M = 0`), polar (parity `(-1)^L`) sector the `+` spinor perturbation is
  `e^{-i omega t} ( a1 Y chi_up + a2 e^{i phi} Y' chi_dn,  i (sigma.n)(b1 Y chi_up + b2 e^{i phi} Y' chi_dn) )`, `Y = P_L(cos theta)`,
  i.e. the `kappa = -(L+1)` and `kappa = L` channels. The `-` spinor follows from the `R_x(pi)` x flavour-swap symmetry (the relative sign of the spin-flip part is fixed by requiring no axial parts in `T_ab`; the other sign mixes in the `L +- 1` axial sectors).
- **Metric.** Regge-Wheeler gauge, `h_tt = alpha^2 H0 Y`, `h_tr = H1 Y`, `h_rr = X^2 H2 Y`, `h_AB = r^2 K Y Omega_AB`, symmetric tetrad perturbation. Unlike for scalar and Proca stars, `H0 != H2`: the traceless angular Einstein equation is sourced by the spin-flip amplitudes `a2`, `b2` (anisotropic stress at linear order), so the single-function ansatz of arXiv:2501.05342 does not apply to Dirac stars.
- **Method.** General `L` through `Y'' = -cot Y' - L(L+1) Y`, with the projections evaluated at generic rational points on the sphere. Output: `relaxation/dirac_nonradial_coeffs.py` (8 Dirac rows and all 7 Einstein components).
- **Checks** (in the script): the background equations, closure of the angular projections at two points, no axial components, the `-` spinor as the symmetry image of the `+` one, and a pure-gauge test (Kosmann Lie derivative; slow, about an hour per sign convention). The `-` sign fails; the run for the `+` sign was stopped before it finished, so this test is still open (the analogous test passes for Proca stars).
- **Independent checks** of the tables (`relaxation/nonradial_checks/`):
  - the projected Bianchi identities hold to finite-difference accuracy in the exterior for arbitrary metric perturbations, and inside the star with on-shell matter;
  - static vacuum solutions are reproduced (flat tidal field; Schwarzschild `L = 2`, `H = r (r - 2M)`, `K = r^2 - 2M^2`);
  - for `L = 0` (with `K = H1 = 0`) the tables reproduce the radial spectrum of `dirac_modes.py` to six digits (e.g. `sigma = 0.017851, 0.034784, 0.042215` at `A_central = 0.05`; unstable 0.02823 vs 0.02824 at 0.25).

`relaxation/dirac_nonradial_modes.py` (`NonRadialModes`, `classified`) solves the eigenproblem; see its docstring for the formulation. Its main choices:
- Einstein components imposed: `Eb, Erth, Etth, Etr`. The Bianchi identities then give `Ett = Ea = 0` for `lambda != 0`. Using the Hamiltonian constraint instead leaves two constraint-violating branches and spurious growing modes.
- Centre: the leading power `r^L` (`r^(L-1)` for `b`, `H1`) is factored out.
- Outer boundary: Robin conditions `F' + (lambda X/alpha + (L+1)/r) F = 0` for the metric. Dirichlet conditions spoil the interior constraints through the spectral derivative.
- An eigenvalue counts as physical only if it converges between two resolutions and satisfies the three un-imposed Einstein components (relative residual < 1e-2).
- `scan_nonradial.py <binary> <out.pkl> <A,...> <L,...> [R] [N1] [N2]` scans single stars of the ground-state family.
- `L = 1` is not reliable with this system (residual gauge freedom). The translation zero mode is recovered with the overdetermined pencil (`eigen_projected(full=False)`).

Results (`G = mu = 1`, `R = 60–100`, `N` up to 190):
- **Ground state: stable for `L = 2, 3, 4`** along the whole family scanned (`A_central` 0.02–0.30, past the maximum mass at 0.23). No converged mode has `Re lambda > 0`.
  - The lowest bound `L = 2` frequency rises from 0.050 (`A_central = 0.08`) to 0.121 (0.20), 0.134 (0.225) and 0.171 (0.30). It stays at about 0.8 of the continuum edge `mu - omega`, so there is no zero mode (no bifurcation to non-spherical solutions).
  - Bound `L = 3` modes: 0.061 (0.08), 0.082 (0.11), 0.1415 (0.20). Bound `L = 4` modes: 0.066 (0.08), 0.162 (0.225).
  - Audits at `A_central = 0.20` (`L = 3`) and 0.225 (`L = 4`) over five `(R, N)` combinations find no growing mode with `|lambda| < 1`.
  - The only "unstable" eigenvalues are artefacts that fail the checks: a high-frequency band near `|Im lambda| = mu + omega`, and continuum box modes with `Re lambda < 1e-3` whose frequencies move with `R`.
- **First excited state: unstable to `L = 2`, with a zero mode at `omega = 0.850`.** The `L = 2` spectrum has a real (non-oscillatory) unstable pair, so these stars also break spherical symmetry.

  | `A_central` | 0.025 | 0.05 | 0.11 | 0.17 | 0.20 | 0.22 | 0.24 | 0.25 | 0.255 |
  |---|---|---|---|---|---|---|---|---|---|
  | `omega` | 0.9779 | 0.9575 | 0.9153 | 0.8821 | 0.8687 | 0.8609 | 0.8538 | 0.8506 | 0.8491 |
  | `L = 2` `lambda` | 0.0024 | 0.0047 | 0.0092 | 0.0111 | 0.0104 | 0.0090 | 0.0061 | 0.0027 | stable (`sigma = 0.0031`) |

  - At small amplitude `lambda` is about `0.095 A_central`, the Newtonian scaling, so the instability extends into the Newtonian limit. A real pair comes with a stable pair of similar magnitude (for example `+-0.0047` and `+-0.0043 i` at 0.05). This is consistent with gravity splitting the near-degenerate 2s and 1d levels of the nearly harmonic core potential.
  - `lambda^2` crosses zero at `A_central = 0.252`, `omega = 0.850`, `M = 1.387`, just past the excited family's maximum mass. This is the Dirac analogue of the `l = 2` bifurcation of one-node scalar boson stars (`omega = 0.858` in arXiv:2501.05342). There the non-spherical (chain/ring) families branch off. Here the spherical stars remain radially unstable on both sides of the crossing.
  - `L = 3, 4` at `A_central = 0.05, 0.14, 0.25`: stable (lowest frequencies 0.0098 / 0.0188, 0.031 / 0.054, 0.061 / 0.094), so `L = 2` is the only non-radial instability found for the one-node stars.

### Non-axisymmetric modes and the axial sector of free Dirac stars

**Why `M = 0` is enough.** The free two-spinor system is invariant under a rotation combined with a U(2) flavour rotation of `(psi_+, psi_-)`: both spinors have the same mass, and the Soler term `(sum psibar psi)^2` would not break this either. The background is invariant under these combined rotations. The linearised operator therefore commutes with them, and its modes form `(2L+1)`-plets labelled by total `(L, M)` and parity. A non-axisymmetric (`M != 0`, e.g. `m = 1`) mode is a rotation of an `M = 0` mode with the same `L`, parity and eigenvalue. The polar `M = 0` analysis above covers parity `(-1)^L`; the axial sector (parity `(-1)^(L+1)`) below completes the non-radial spectrum. (The `L = 0` axial sector does not couple to the metric, so its modes are those of the Dirac operator on the fixed background and are stable.)

`derivations/dirac_axial.py` derives the axial `M = 0` equations with the machinery of `dirac_nonradial.py`:
- **Ansatz.** `sigma.n` maps the polar channels `kappa = -(L+1), L` onto `kappa = L+1, -L` (the same `j'`, opposite parity). The perturbation is therefore the polar one with upper and lower components exchanged: `dpsi_+ = e^{-i omega t}((sigma.n)(c1 Y chi_up + c2 e^{i phi} Y' chi_dn), i (d1 Y chi_up + d2 e^{i phi} Y' chi_dn))`.
- **Metric.** RW axial, `h_tphi = h0 sin Y'`, `h_rphi = h1 sin Y'`.
- **Flavour sign.** The relative sign of the `-` spinor's spin-flip part is again `s2 = -1`; `s2 = +1` leaves polar components in `T_ab`.
- **Checks.** The projections close, the polar Einstein components vanish, and the `-` spinor is the symmetry image of the `+` one. Only the spin-flip channel `(c2, d2)` sources the metric.
- **Bianchi identity.** The `phi` component, `-(lambda Etph - (alpha alpha'/X^2) Erph)/alpha^2 + (Erph' - (X'/X + 1/r) Erph)/X^2 + 3 Erph/(r X^2) - (L(L+1) - 2) Ethph/r^2 = 0`, holds to 1e-8 for `L = 1, 2, 3` in the Schwarzschild exterior (arbitrary `h0`, `h1`) and inside a star with on-shell matter (`relaxation/nonradial_checks/bianchi_ax.py`).

`relaxation/dirac_axial_modes.py` (`AxialModes`, `classified`):
- **Equations imposed.** The 8 Dirac rows and `Etph, Erph` (`L >= 2`; the Bianchi identity then forces `Ethph = 0`, which is used as the check), or `Etph` with the gauge `h1 = 0` at `L = 1` (`Erph` is the check).
- **Discretisation.** Regular factor `p = L - 1, L, L + 1` for `c`, `d`, `h0/h1`; Robin outer condition `F' + (lambda X/alpha + L/r) F = 0` for the metric.
- **Zero mode.** At `L = 1` the solver recovers the exact zero mode predicted by the symmetry, the relative phase rotation `dpsi_+ = i psi_+`, `dpsi_- = -i psi_-` (`|lambda| < 1e-5`; the eigenvector is purely in the phases of both channels, with no metric perturbation at leading order).

**Polar `L = 1` revisited.** At `L = 1`, `Y'' = -Y`, so only the combination `Ea - Eb` of the angular components exists. Checking `Ea` and `Eb` separately (as before) wrongly flagged every `L = 1` mode; `dirac_nonradial_modes.residuals` now uses `Ea - Eb`. With this, the converged `L = 1` modes are stable:
- ground state: lowest `sigma` 0.0198, 0.0424, 0.0532, 0.0636, 0.0737, 0.1263 at `A_central` = 0.05, 0.11, 0.14, 0.17, 0.20, 0.225;
- excited state: lowest `sigma` 0.0057, 0.0108, 0.0125, 0.0205 at `A_central` = 0.05, 0.11, 0.17, 0.25.

The only growing candidate is a real eigenvalue that drifts with resolution and domain (0.048, 0.069, 0.090 at `N` = 110, 150, 190 for `A_central = 0.05`) and violates the constraints at the 0.2–0.5 level. It is presumably the residual `L = 1` gauge freedom (accelerated translations) and is not physical.

**Results.** The axial sector is stable for `L = 1, 2, 3` along the ground-state family (`A_central` 0.02–0.30) and the first-excited family (0.025–0.25):
- Apart from the flavour zero mode, no converged mode has `Re lambda > 5e-4`.
- Lowest axial frequencies: ground state `L = 1`: 0.027, 0.042, 0.070, 0.082, 0.094, 0.104, 0.113, 0.129 (`A_central` = 0.05 ... 0.30); `L = 2`: 0.020, 0.033, 0.057, 0.069, 0.081, 0.091, 0.101 (0.05 ... 0.25).
- Excited state, near-zero but stable: `L = 1` 0.0020 (0.025), 0.0041 (0.05); `L = 2` 0.0029, 0.0055; `L = 3` 0.0021, 0.0045.
- The only rejected candidates are the high-frequency (`|Im lambda| ~ 0.6–0.9`) unconverged artefacts.

Hence, for free Dirac stars, **there is no `m = 1` (or any `M != 0`) instability**: the non-radial instabilities are exactly those of the polar `M = 0` analysis (none for the ground state; `L = 2` for the excited state, with its `M = +-1, +-2` copies). For Proca stars the same symmetry argument (plain rotations) implies that the polar `L = 1` oscillatory instability found below also has `M = +-1` copies, i.e. an `m = 1` non-axisymmetric instability with the same eigenvalue.

### Non-radial (axisymmetric) perturbations of free Proca stars

`derivations/proca_nonradial.py` derives the linearised polar (`M = 0`, parity `(-1)^L`) equations for free Proca stars, in ProcaStar's conventions (`A = e^{-i omega t}(f dt + i g dr)`):
- **Perturbation.** `dA_t = e^{-i omega t} a0 Y`, `dA_r = e^{-i omega t} a1 Y`, `dA_theta = e^{-i omega t} a2 Y'`. The metric is in the same Regge-Wheeler form as for Dirac stars.
- **Anisotropic stress.** The traceless angular Einstein equation is not sourced by the matter at linear order (`F_theta r` and `A_theta` vanish in the background), so `H0 = H2`, as assumed in arXiv:2501.05342.
- **Checks** (in the script): the stress tensor reproduces ProcaStar's `rho`; the background equations hold; the projections close, with no `phi` or axial components; and the pure-gauge test (`h = L_zeta g`, `dA = L_zeta A`, general polar `zeta`) is satisfied identically.
- **Coefficients** are written to `relaxation/proca_nonradial_coeffs.py`.

`relaxation/proca_relax.py` re-solves ProcaStar's shooting solutions by relaxation (free field). The family has `M_max = 1.0581` at `A_central = 0.0375`, `omega = 0.874`.

`relaxation/proca_nonradial_modes.py` (`ProcaModes`, `classified`) uses the formulation of the Dirac solver: Einstein components `Eb, Erth, Etth, Etr`, regularity factored out at the centre, Robin outer conditions, and classification by convergence plus the residuals of the un-imposed components.
- The `lambda^2` terms of the `r` and `theta` Proca equations are linearised with auxiliary fields.
- At `L = 1`, `K = 0` and `Eb` is dropped; the angular check uses `Ea - Eb`, the only combination present when `Y'' = -Y`.
- **L = 0 check.** With `K = H1 = 0` the fundamental radial frequency converges (`sigma_1` = 0.01408, 0.01495, 0.01005, 0.00582 at `A_central` = 0.01, 0.025, 0.033, 0.036) and reaches zero between 0.036 and 0.039, at the maximum mass.
- Drivers and checks are in `relaxation/proca_nonradial_checks/`.

Results (`G = mu = 1`):
- **L = 2: real (non-oscillatory) instability on the whole stable branch, from the Newtonian limit to just past the maximum mass.** The zero mode is at `omega = 0.8713` (`A_central = 0.0391`, `M = 1.058`). This is the `l = 2` bifurcation found for Proca stars in arXiv:2501.05342 (`omega = 0.8713`, `M = 1.0579`) with the static zero-mode method. Beyond it the `L = 2` mode is stable (`sigma` = 0.0015, 0.0022, 0.0062, 0.0090, 0.0139 at `A_central` = 0.0395, 0.040, 0.045, 0.05, 0.06). This is the quadrupolar deformation towards the prolate stars of arXiv:2311.14800. The fastest e-folding time is about 160 / mu, near `omega = 0.89`.

  | `A_central` | 0.003 | 0.006 | 0.01 | 0.015 | 0.02 | 0.025 | 0.03 | 0.033 | 0.036 | 0.038 | 0.039 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | `omega` | 0.9671 | 0.9502 | 0.9334 | 0.9173 | 0.9044 | 0.8938 | 0.8847 | 0.8799 | 0.8755 | 0.8727 | 0.8714 |
  | `L = 2` `lambda` | 0.00174 | 0.00309 | 0.00448 | 0.00566 | 0.00629 | 0.00638 | 0.00584 | 0.00511 | 0.00387 | 0.00239 | 0.00069 |

  At small amplitude `lambda` is about 0.06 `(mu - omega)`, the Newtonian scaling. The values are converged in `N` (to 190) and independent of `R` (0.00629 at `R` = 60 and 90 for `A_central = 0.02`), with residuals 1e-4 to 1e-3.
- **L = 1: oscillatory instability from `omega` of about 0.928 onwards**, continuing past the maximum mass. It sets in through a Krein collision: a single stable pair at `sigma = 0.0197` for `A_central = 0.011` becomes the quartet `+-0.0008 +- 0.0204 i` at 0.012.

  | `A_central` | 0.012 | 0.015 | 0.02 | 0.025 | 0.03 | 0.036 | 0.04 | 0.05 |
  |---|---|---|---|---|---|---|---|---|
  | `L = 1` `lambda` | 0.0008 + 0.0204 i | 0.0018 + 0.0230 i | 0.0028 + 0.0265 i | 0.0037 + 0.0295 i | 0.0043 + 0.0321 i | 0.0050 + 0.0348 i | 0.0054 + 0.0363 i | 0.0060 + 0.0397 i |

  - The eigenvalue agrees to five digits between `N` = 110 and 150 and between `R` = 60 and 90. The Hamiltonian constraint holds to 1e-7 or better in the bulk.
  - Caveat: in the RW-type gauge with `K = 0`, `L = 1` keeps a residual gauge freedom (time-dependent translations). The eigenvector carries a component of it (`H0` grows roughly linearly outwards). A gauge-invariant reformulation or an evolution should confirm this mode.
  - In the Newtonian picture the spherical star is a "hedgehog" `A_i ~ g(r) n_i`. The perturbation `dA ~ epsilon_z psi(r)` towards a fixed-polarisation s-wave state is exactly polar `L = 1`. An oscillatory instability in this channel is the analogue of the decay of an excited state.
- **L = 3: stable** along the family (lowest frequencies 0.014, 0.027, 0.039, 0.048, 0.053, 0.060 at `A_central` = 0.006, 0.015, 0.025, 0.033, 0.038, 0.045).

The momentum-constraint diagnostic `Mom_Z` had a sign error in its `A_ww h_ww'` term, now fixed. In dynamical runs it produced a resolution-independent residual. It now converges at 4th order and agrees with the constraint computed independently from the slices. Values of `Mom_L2` from earlier runs should not be relied on.

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
