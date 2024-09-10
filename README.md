C++ code that constructs and evolves spherical boson star solutions in polar-areal gauge, optionally converts to isotropic gauge, then evolves using the BSSN formalism reduced to spherical symmetry using a modified cartoon approach.
Can also generate families of BS models if the cycle_models line is uncommented in main(). More detailed rundown below:

NOTE: for some reason  my .gitignore didn't seem to work. Nothing ending in .dat or .cbp in the repo is needed to run the code (Can look at the .gitignore for a full list)

###COMPILE + RUN####
I include a shell script make_BSEvolver_gcc.sh that should produce the executable BosonStarStability.o in the directory it is run using gcc (the shell script must be run in the same directory as all source files currently).
The executable needs to run in the same directory as a BSParams.par file (included in repo), else it will abort. Many of these are self-explanatory/same as GRChombo but some I will explain below.

###OUTPUTS###
Outputs all end in .dat. The files that are produced (in order of a usual run):
-BSdata.dat: stores field values from original BS construction, in order: r, A, X, phi, t

-isotropic.dat: if "isotropic" is true in params file, spits out isotropic field values; this is messy + probably won't need but can check source code write_isotropic() for order

-SliceData.dat: written to every write_interval timesteps. Stores z-value followed by the 13 BSSN + mater variables, whose order is given in the enum class bssn_var at the top of EvolutionVariables.h as well as in the struct below it.
I try to be consistent with the ordering of these variables throughout.

-Diagnostics.dat: stores z then profiles of  Ham + Momentum constraints, determinant of h (conformal metric), auxiliary constraint, then some test values that I set in source code

-constraint_norms.dat: written to every write_CN_interval timesteps during evolve(); stores L^2 norm of Ham + Mom over time, then central values of chi + field amplitude, then central energy density.

There is also conv.dat only written to as the result of a convergence test (currently only fully implemented on the BS solver) and BosonStars.dat written to to store the results of model cycling, both commented out in main by default.

###BASIC OVERVIEW OF STRUCTURE###
The above should be enough to use the code in principle, but here's a brief overview of what it's doing under the hood...

First, a BosonStar object is created. This is associated with an array of the variables used in areal-radius gauge model construction, as well as grid/potential parameters that are needed during initial data construction.
There are also arrays for the isotropic values that will be separately allocated and filled if isotropic is true. All of this is governed by the BosonStars header + .cpp file.

This is then translated into something we can dynamically evolve. The structure of the classes in EvolutionVariables.h is as follows. At the bottom level is the BSSNState struct, an instance of which contains a set of 
BSSN (+ matter) variables at a point. This is useful as we can e.g. overload arithmetic operators, so adding two BSSNStates adds the variables elementwise, saving a lot of rewriting.
Next is the BSSNSlice class, which holds an array of BSSNStates corresponding to a time slice and some slice-relevant info. Partial derivatives are defined here, so we can take them on a slice without associating it with a full spacetime.
Finally there is the spacetime class, which holds an array of slices and a bunch of arrays holding auxilary quantities used in the BSSN equations so they do not need to be re-computed, but which get overwritten on every time step.
The main functions relevant here are initialize, which takes a BosonStar object and fixes the initial time slice, and evolve, which runs dynamical evolution-- anything else could be made private if i didn't want it accessible for frequent
testing in debugging. Initialize must be run first or UB may happen. The integer max_stored_slices is essentially the maximum size of the slices array in a spacetime; once past it the array will cycle back by 1, overwriting the smallest
element on each timestep.

Another important thing to note is the current_slice_ptr, which may be unintuitive. Since I want partial derivatives accessible on a slice but also to be callable as members of the spacetime, each spacetime has a pointer
current_slice_ptr which determines on which slice derivatives are evaluated (unless overruled by an optional argument to d_z). This must therefore be manually updated to hold the address of the desired slice 
before anything evolution-related can be done; see e.g. the RK4 implementation in evolve() for example. 

