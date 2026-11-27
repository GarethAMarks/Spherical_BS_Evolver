# MPI build: requires an MPI C++ compiler wrapper (mpicxx), e.g. from Open MPI or MPICH.
# Produces BosonStarStability_mpi.o; run with e.g. mpirun -np 4 ./BosonStarStability_mpi.o
mkdir -p obj
mkdir -p obj/MPI/finite_temp
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "mathutils.cpp" -o obj/MPI/mathutils.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "BosonStar.cpp" -o obj/MPI/BosonStar.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "LinearPerturbation.cpp" -o obj/MPI/LinearPerturbation.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "EvolutionVariables.cpp" -o obj/MPI/EvolutionVariables.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "main.cpp" -o obj/MPI/main.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "ComplexScalarField.cpp" -o obj/MPI/ComplexScalarField.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "RealScalarField.cpp" -o obj/MPI/RealScalarField.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "FiniteTempBackground.cpp" -o obj/MPI/FiniteTempBackground.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "FiniteTempMatter.cpp" -o obj/MPI/FiniteTempMatter.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "finite_temp/FiniteTempRelaxer.cpp" -o obj/MPI/finite_temp/FiniteTempRelaxer.o
mpicxx -Wall -fexceptions -g -DUSE_MPI -I.  -c "finite_temp/seed_shooter.cpp" -o obj/MPI/finite_temp/seed_shooter.o
mpicxx  -o BosonStarStability_mpi.o obj/MPI/FiniteTempBackground.o obj/MPI/FiniteTempMatter.o obj/MPI/finite_temp/FiniteTempRelaxer.o obj/MPI/finite_temp/seed_shooter.o obj/MPI/BosonStar.o obj/MPI/LinearPerturbation.o obj/MPI/EvolutionVariables.o obj/MPI/main.o obj/MPI/mathutils.o obj/MPI/ComplexScalarField.o obj/MPI/RealScalarField.o  -lm
