mkdir obj
mkdir obj/Debug
g++ -Wall -fexceptions -g -I.  -c "mathutils.cpp" -o obj/Debug/mathutils.o
g++ -Wall -fexceptions -g -I.  -c "BosonStar.cpp" -o obj/Debug/BosonStar.o
g++ -Wall -fexceptions -g -I.  -c "LinearPerturbation.cpp" -o obj/Debug/LinearPerturbation.o
g++ -Wall -fexceptions -g -I.  -c "EvolutionVariables.cpp" -o obj/Debug/EvolutionVariables.o
g++ -Wall -fexceptions -g -I.  -c "main.cpp" -o obj/Debug/main.o
g++ -Wall -fexceptions -g -I.  -c "ComplexScalarField.cpp" -o obj/Debug/ComplexScalarField.o
g++  -o BosonStarStability.o obj/Debug/BosonStar.o obj/Debug/LinearPerturbation.o obj/Debug/EvolutionVariables.o obj/Debug/main.o obj/Debug/mathutils.o obj/Debug/ComplexScalarField.o   -lm
