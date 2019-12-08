#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include <mpi.h>
#include "quantum_dots.hpp"


void grid_search_psi1(Psi1 instance, double start, double nudge, int steps)
{
  double min_en  = 1e6;
  double min_par = 0;
  instance.reset();
  instance.varparam = start;

  for(int i = 0; i<steps; i++)
  {
    instance.metropolis();

    if(instance.averages(instance.mc_cycles, 0) < min_en)
    {
      min_en = instance.averages(instance.mc_cycles, 0);
      min_par = instance.varparam;
    }
    instance.reset();
    // reset average and positions
    instance.varparam += nudge; // change variational parameter slightly
  }
  std::cout << min_par << " En: "<< min_en << std::endl;
}

void grid_search_psi2(Psi2 instance, double start, double nudge, int steps)
{
  double min_en  = 1e6;
  double min_par = 0;
  instance.reset();
  instance.varparam = start;

  for(int i = 0; i<steps; i++)
  {
    instance.metropolis();

    if(instance.averages(instance.mc_cycles, 0) < min_en)
    {
      min_en = instance.averages(instance.mc_cycles, 0);
      min_par = instance.varparam;
    }
    instance.reset();
    // reset average and positions
    instance.varparam += nudge; // change variational parameter slightly
  }
  std::cout << min_par << " En: "<< min_en << std::endl;
}

int main(int argc, char *argv[])
{
  int num_processes;
  int process_rank;
  int n = 1000000;
  double time_start = 0;
  double execution_time = 0;
  double freq = 1;
  double alpha_init = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  time_start = MPI_Wtime();
  Psi1 model(n, alpha_init, freq, process_rank); // n, alpha, omega, seed

  model.metropolis();
  std::cout << model.averages(100000, 0) << std::endl;
  execution_time = MPI_Wtime() - time_start;

  grid_search_psi1(model, 0.1 + process_rank, 0.1, 50);
  MPI_Finalize();
  return 0;
}

// TODO: Overload QD constructor to accept both psi1 and psi2 objects!! Check acceptance rate
