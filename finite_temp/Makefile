# Standalone build of the finite-temperature relaxer (as in the original Finite_Temp_Relaxer repo).
# SBSE itself compiles these sources via ../make_BSEvolver_gcc.sh; this Makefile is only for running the
# relaxer on its own, e.g. to produce profile.dat for SBSE's finite_temp = 1, ft_relax = 0 mode.
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2
CPPFLAGS ?= -I..
LDLIBS ?= -lm

TARGET := finite_temp_relaxer
SOURCES := relaxer_main.cpp FiniteTempRelaxer.cpp ../mathutils.cpp seed_shooter.cpp
OBJECTS := relaxer_main.o FiniteTempRelaxer.o mathutils.o seed_shooter.o
DEPFILES := $(OBJECTS:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

mathutils.o: ../mathutils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPFILES)

clean:
	rm -f $(OBJECTS) $(DEPFILES) $(TARGET)
