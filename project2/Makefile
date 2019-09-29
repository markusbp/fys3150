CC = g++
CXXFLAGS := -Wall -g -o
LIB = -larmadillo
all: conduct_unit_tests main_file

conduct_unit_tests: unit_tests.cpp
	$(CC) $(CXXFLAGS) conduct_unit_tests unit_tests.cpp tools.cpp jacobi_rotation.cpp $(LIB)

main_file: main.cpp
			$(CC) $(CXXFLAGS) main_file main.cpp tools.cpp jacobi_rotation.cpp $(LIB)
