# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC = g++
CXXflags = -Wall -o
# Here we define the library functions we nee
LIB = -larmadillo

all: analysis parallel_ising run_unit_tests

analysis: data_analysis.o ising_model.o
	 		$(CC) $(CXXflags) analysis data_analysis.o ising_model.o $(LIB)

parallel_ising: parallel_ising.cpp ising_model.cpp
			mpicxx -O2 -o parallel_ising ising_model.cpp parallel_ising.cpp $(LIB)

run_unit_tests: unit_tests.o ising_model.o
			$(CC) $(CXXflags) run_unit_tests unit_tests.o ising_model.o  $(LIB)

data_analysis.o: data_analysis.cpp
			$(CC) -c data_analysis.cpp $(LIB)

ising_model.o: ising_model.cpp
			$(CC) -c ising_model.cpp $(LIB)

run_unit_tests.o: unit_tests.cpp
			$(CC) -c unit_tests.cpp $(LIB)

clean:
	rm -f core ising_model.o data_analysis.o run_unit_tests.o

#run_unit_tests: unit_tests.o
#			$(CC) $(CXXflags) unit_tests.cpp
