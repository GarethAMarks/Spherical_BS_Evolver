g++ -Wall -fexceptions -g  -c "mathutils.cpp" -o obj/Debug/mathutils.o
g++ -Wall -fexceptions -g  -c "BosonStar.cpp" -o obj/Debug/BosonStar.o
g++ -Wall -fexceptions -g  -c "EvolutionVariables.cpp" -o obj/Debug/EvolutionVariables.o
g++ -Wall -fexceptions -g  -c "/home/gmarks2000/Documents/codeBlocks projects/BosonStarStability/main.cpp" -o obj/Debug/main.o
g++  -o bin/Debug/BosonStarStability obj/Debug/BosonStar.o obj/Debug/EvolutionVariables.o obj/Debug/main.o obj/Debug/mathutils.o   -lm
