# MPI build: requires an MPI C++ compiler wrapper (mpicxx), e.g. from Open MPI or MPICH.
# Produces BosonStarStability_mpi.o; run with e.g. mpirun -np 4 ./BosonStarStability_mpi.o
mkdir -p obj
mkdir -p obj/MPI
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "mathutils.cpp" -o obj/MPI/mathutils.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "BosonStar.cpp" -o obj/MPI/BosonStar.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "LinearPerturbation.cpp" -o obj/MPI/LinearPerturbation.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "EvolutionVariables.cpp" -o obj/MPI/EvolutionVariables.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "main.cpp" -o obj/MPI/main.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "ComplexScalarField.cpp" -o obj/MPI/ComplexScalarField.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "RealScalarField.cpp" -o obj/MPI/RealScalarField.o
mpicxx  -o BosonStarStability_mpi.o obj/MPI/BosonStar.o obj/MPI/LinearPerturbation.o obj/MPI/EvolutionVariables.o obj/MPI/main.o obj/MPI/mathutils.o obj/MPI/ComplexScalarField.o obj/MPI/RealScalarField.o  -lm
