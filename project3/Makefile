# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC = g++
INC_DIR = ../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary
CXXflags = -Wall -g -I$(INC_DIR)
# Here we define the library functions we nee
LIB = -larmadillo
# Here we define the name of the executable
all :  gauss_quadrature monte_carlo parallel_monte_carlo

#unit_tests :	   	unit_tests.o
#			${CXXFLAGS} unit_tests.o ${LIB} -o unit_tests

gauss_quadrature : gauss_quadrature.o lib.o gauss_laguerre.o
		 $(CC)	${CXXFLAGS} gauss_quadrature.o lib.o gauss_laguerre.o ${LIB} -o gauss_quadrature

monte_carlo : monte_carlo.o
			$(CC) ${CXXFLAGS} monte_carlo.o ${LIB} -o monte_carlo

parallel_monte_carlo : parallel_monte_carlo.cpp
			mpicxx -O2 -o parallel_monte_carlo parallel_monte_carlo.cpp -larmadillo

#unit_tests.o :			unit_tests.cpp
#		        	${CXXFLAGS} -c unit_tests.cpp

lib.o : $(INC_DIR)/lib.cpp
			$(CC) $(CXXFLAGS) -c $(INC_DIR)/lib.cpp

gauss_laguerre.o : $(INC_DIR)/gauss_laguerre.cpp
			$(CC) $(CXXFLAGS) -c $(INC_DIR)/gauss_laguerre.cpp

gauss_quadrature.o :	gauss_quadrature.cpp
		        $(CC) ${CXXFLAGS} -c gauss_quadrature.cpp  $(INC_DIR)/gauss_laguerre.cpp $(INC_DIR)/lib.cpp

monte_carlo.o : monte_carlo.cpp
						$(CC)	${CXXFLAGS} -c monte_carlo.cpp

clean:
	 		rm -f core monte_carlo.o gauss_quadrature.o lib.o gauss_laguerre.o unit_tests.o
