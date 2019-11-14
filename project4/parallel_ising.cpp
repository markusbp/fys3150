#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include <mpi.h>
#include "ising_model.hpp"


int main(int argc, char *argv[])
{
  int l = atoi(argv[1]); // Size of spin grid, lxl
  int n = atoi(argv[2]); // number of Monte Carlo cycles (each of lxl spin flips)
  double t; // temperature
  double dt;
  int num_processes = 0;
  int process_rank = 0;

  // init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD); // Distribute l,n, t to all procs
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Do parallel calculation with many temperatures! Remember timing

  dt = 0.05;

  t = 2.0 + dt*process_rank; // 2.0 2.05, 2.1, 2.15 ... 2.020 2.025, 2.030 2.035

  IsingModel model(l, n, t);

  model.set_random_spin_config(process_rank); //think about seed!
  model.metropolis(process_rank);
  std::cout << model.mean_energy(n-1) << std::endl;
  MPI_Finalize();


}
