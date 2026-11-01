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
- bianchi_ax.py    phi component of the Bianchi identity for the axial tables (vacuum and on-shell interior)
- ax1.py, ax0.py   axial spectra of single stars (python3 ax1.py A L R [exc]); inspection of near-zero eigenvalues
- l1d.py, l1aud.py polar L = 1 spectra along a family; resolution / domain audit of L = 1 candidates
