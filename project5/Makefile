# Here we define compiler option, libraries and the  target
CC = g++
CXXflags = -Wall -o
# Here we define the library functions we nee
LIB = -larmadillo

all: run_simulation run_unit_tests search_psi1

run_simulation: simulate_qd.cpp quantum_dots.cpp
		mpicxx -O2 -o run_simulation simulate_qd.cpp quantum_dots.cpp $(LIB)

run_unit_tests: unit_tests.o
		$(CC) $(CXXflags) run_unit_tests unit_tests.o $(LIB)

search_psi1: psi1.o quantum_dots.o
		$(CC) $(CXXflags) search_psi1 psi1.o quantum_dots.o $(LIB)

quantum_dots.o: quantum_dots.cpp
		$(CC) -c -O2 quantum_dots.cpp $(LIB)

psi1.o: psi1.cpp quantum_dots.cpp
		$(CC) -c -O2 psi1.cpp quantum_dots.cpp $(LIB)

unit_tests.o: unit_tests.cpp
		$(CC) -c -O2 unit_tests.cpp  $(LIB)

clean:
	rm -f core unit_tests.o
