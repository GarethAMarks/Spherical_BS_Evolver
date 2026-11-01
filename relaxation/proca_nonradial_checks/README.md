Drivers and checks for the non-radial Proca-star modes (main README, "Non-radial (axisymmetric) perturbations of free
Proca stars"). SBSE_BINARY: the compiled evolver (default ../../BosonStarStability); SBSE_WORK: scratch directory for
cached backgrounds (default .).
- fam.py       backgrounds (ProcaStar shooting + proca_relax), cached; prints the family
- l0.py        L = 0 reduction: fundamental radial frequency -> 0 at the maximum mass
- t1.py, scan.py   classified spectra for one star / along the family (python3 scan.py A1,A2 L N1,N2 out.pkl [R])
- audit1.py, checks.py, ettprof.py   resolution / domain audits, small-amplitude and crossing checks, constraint profile
- figp.py      summary figure
