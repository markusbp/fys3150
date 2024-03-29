# Here we define compiler option, libraries and the  target
CC = g++
CXXflags = -Wall -o
# Here we define the library functions we nee
LIB = -larmadillo

all: run_unit_tests virial search_psi1 search_psi2 run_psi2

search_psi2: search_psi2.cpp quantum_dots.cpp
		mpicxx -O2 -o search_psi2 search_psi2.cpp quantum_dots.cpp $(LIB)

run_unit_tests: unit_tests.o quantum_dots.o
		$(CC) $(CXXflags) run_unit_tests unit_tests.o quantum_dots.o $(LIB)

virial: virial.o quantum_dots.o
		$(CC) $(CXXflags) virial virial.o quantum_dots.o $(LIB)

search_psi1: psi1.o quantum_dots.o
		$(CC) $(CXXflags) search_psi1 psi1.o quantum_dots.o $(LIB)

run_psi2: psi2.o quantum_dots.o
		$(CC) $(CXXflags) run_psi2 psi2.o quantum_dots.o $(LIB)

quantum_dots.o: quantum_dots.cpp
		$(CC) -c -O2 quantum_dots.cpp $(LIB)

psi1.o: psi1.cpp quantum_dots.cpp
		$(CC) -c -O2 psi1.cpp quantum_dots.cpp $(LIB)

psi2.o: psi2.cpp quantum_dots.cpp
		$(CC) -c -O2 psi2.cpp quantum_dots.cpp $(LIB)

virial.o: virial.cpp quantum_dots.cpp
		$(CC) -c -O2 virial.cpp quantum_dots.cpp $(LIB)

unit_tests.o: unit_tests.cpp
		$(CC) -c -O2 unit_tests.cpp  $(LIB)

clean:
	rm -f core unit_tests.o psi1.o search_psi1.o psi2.o search_psi2.o quantum_dots.o virial.o
