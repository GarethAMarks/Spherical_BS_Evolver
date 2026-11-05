# Serial build with optimisation. Produces BosonStarStability.o; run with ./BosonStarStability.o next to BSParams.par.
# Object files go to obj/Release, so this build can coexist with the MPI build (obj/MPI).
# Override the flags with e.g.  OPTFLAGS="-O0 -g" bash make_BSEvolver_gcc.sh  (debugging)
#                          or  OPTFLAGS="-O3 -march=native" bash make_BSEvolver_gcc.sh  (machine-specific binary).
# Do not add -ffast-math: the solvers rely on NaN checks (isnan) that it disables.
set -e

CXX=${CXX:-g++}
OPTFLAGS=${OPTFLAGS:-"-O2"}
CXXFLAGS="-Wall -fexceptions $OPTFLAGS -I."
OBJDIR=obj/Release

SOURCES="mathutils BosonStar LinearPerturbation EvolutionVariables main ComplexScalarField RealScalarField"

mkdir -p "$OBJDIR"
OBJECTS=""
for src in $SOURCES; do
    echo "compiling $src.cpp"
    $CXX $CXXFLAGS -c "$src.cpp" -o "$OBJDIR/$src.o"
    OBJECTS="$OBJECTS $OBJDIR/$src.o"
done

$CXX $OPTFLAGS -o BosonStarStability.o $OBJECTS -lm
echo "built BosonStarStability.o ($CXX $OPTFLAGS)"
