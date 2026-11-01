Validation and driver scripts for the non-radial Dirac-star modes (see the main README, "Non-radial (axisymmetric)
perturbations of free Dirac stars"). Set SBSE_BINARY to the compiled evolver (default ../../BosonStarStability) and
SBSE_WORK to a scratch directory for shooting runs and cached backgrounds (default .).

- vac.py           static vacuum solutions (flat tidal field, Schwarzschild L = 2) satisfy all 7 Einstein rows
- bianchi_vac.py   projected Bianchi identities, Schwarzschild exterior, arbitrary time-dependent metric perturbation
- bianchi_in.py    the same inside a star, with the matter made on-shell by integrating the 8 Dirac rows
- l0check.py       L = 0 reduction of the tables vs the radial ModeSolver
- bg.py, bg_exc.py ground-state / excited (eigen = 1) backgrounds: shooting + relaxation, cached
- track0.py, pin0.py, audit.py   spectra near zero along the excited family; resolution/domain audits
- fig.py           summary figure
