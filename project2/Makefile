CC = g++
CXXFLAGS := -O3 -o#-Wall -g -o -O3
LIB = -larmadillo
all: conduct_unit_tests harmonic_oscillator

conduct_unit_tests: unit_tests.cpp
	$(CC) $(CXXFLAGS) conduct_unit_tests unit_tests.cpp tools.cpp jacobi_rotation.cpp $(LIB)

harmonic_oscillator: harmonic_oscillator.cpp
			$(CC) $(CXXFLAGS) harmonic_oscillator harmonic_oscillator.cpp tools.cpp jacobi_rotation.cpp $(LIB)
