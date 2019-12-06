#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include <mpi.h>
#include "quantum_dots.hpp"

int main(int argc, char *argv[])
{
  int num_processes;
  int process_rank;
  double time_start = 0;
  double execution_time = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  time_start = MPI_Wtime();
  QuantumDots model(100000, 0.88, 1, 1, process_rank); // n, alpha, beta, omega, seed

  model.metropolis();

  execution_time = MPI_Wtime() - time_start;
  MPI_Finalize();
  std::cout<< model.averages.row(100000-1) << std::endl;
  return 0;
}

// TODO: Overload QD constructor to accept both psi1 and psi2 objects!! Check acceptance rate
