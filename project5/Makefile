# Here we define compiler option, libraries and the  target
CC = g++
CXXflags = -Wall -o
# Here we define the library functions we nee
LIB = -larmadillo

all: run_simulation run_unit_tests

run_simulation: simulate_qd.cpp quantum_dots.cpp
		mpicxx -O2 -o run_simulation simulate_qd.cpp quantum_dots.cpp $(LIB)

run_unit_tests: unit_tests.o
		$(CC) $(CXXflags) run_unit_tests unit_tests.o $(LIB)

unit_tests.o: unit_tests.cpp
		$(CC) -c -O2 unit_tests.cpp $(LIB)

clean:
	rm -f core unit_tests.o