#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include <mpi.h>
#include "quantum_dots.hpp"

// script for performing 2D gridsearch for psi2 under both variational parameters,
// alpha and beta, using MPI to run multiple simulations at once
// usage (example):
// mpirun -np 8 /search_psi2 grid_length grid_step alpha_init beta_init
// grid_length: number of steps in one row of grid search: total steps is grid_length^2
// grid_step: step in parameter, resolution of search
// alpha_init: starting value of alpha
// beta_init: starting value of beta
// -np signifies number of parallel processes. (standard MPI)

int main(int argc, char *argv[])
{
  int process_rank = 0; // for MPI to fill in, gives name/number of process
  int num_processes = 0;
  double time_start = 0;
  double execution_time = 0; // timing

  int n = 2000000; // number of MC cycles, just set to this as default
  int col = 0; // counter in columns to ensure each process does unique work
  int grid_length = atoi(argv[1]); // see top
  double grid_step = atof(argv[2]);
  double alpha_init = atof(argv[3]);
  double beta_init = atof(argv[4]);
  double freq = 1; // harmonic oscillator frequency
  arma::mat variational_params = arma::zeros(2, grid_length); // alpha and beta

  // Initialize arrays/matrices to work with MPI reduce
  double **all_energies;
  double **trial_energies;
  all_energies =  new double*[grid_length];
  trial_energies =  new double*[grid_length];

  for (int i = 0; i < grid_length; i++)
  {
    all_energies[i] = new double[grid_length];
    trial_energies[i] = new double[grid_length];
    variational_params(0, i) = alpha_init + (i+1)*grid_step;
    variational_params(1, i) = beta_init + (i+1)*grid_step; // set all alpha, betas
    for (int j = 0; j < grid_length; j++)
    {
      trial_energies[i][j] = 0; // initialize all elements to zero
      all_energies[i][j] = 0; // initialize all elements to zero
    }
  }

  // init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  time_start = MPI_Wtime();
  // instantiate psi2 quantum dot model, one for each process
  Psi2 trial2(n, alpha_init, beta_init, freq, process_rank);
  // perform parallel grid search
  for(int i = 0; i < grid_length; i++)
  {
    trial2.a = variational_params(0, i); // set alpha
    trial2.varparam = beta_init; // reset beta for each row
    for(int j = 0; j < grid_length; j+= num_processes)
    {
      col = j + process_rank; // gives each process unique beta
      if(i < grid_length && col < grid_length) // don't overstep
      {
        trial2.varparam = variational_params(1, col); // each process gets a beta
        trial2.metropolis(); // run sim
        trial_energies[i][col] = trial2.averages(n, 0); // save energy to find min.
      }
    }
    if(process_rank==0)
    { // keep track of progress
      std::cout << "Computing row " << i+1 << "/" << grid_length << std::endl;
    }
  }

  for (int i = 0; i < grid_length; i++)
  { // Each process has unique elements in data matrix, sum just merges them
    MPI_Reduce(trial_energies[i], all_energies[i], grid_length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  execution_time = MPI_Wtime() - time_start; // I don't want to know :)
  MPI_Finalize();

  if(process_rank == 0)
  {
    // Write result array to file
    std::fstream filewriter;
    std::string name = "./results/psi2/full_gridsearch_psi2";
    filewriter.open(name, std::ios::out); // write mode
    for(int i = 0; i < grid_length; i++)
    {
      for(int j = 0; j < grid_length; j++)
      {
        filewriter << all_energies[i][j];
        if(j < grid_length -1)
        {
          filewriter << ","; // csv file :)
        }
      }
      filewriter << "\n";
    }
    variational_params.save("results/psi2/full_gs_psi2_varparams", arma::csv_ascii);
    filewriter.close();
  }

  for (int i = 0; i < 4; ++i)
  {   // clear memory
      delete [] trial_energies[i];
      delete [] all_energies[i];
  }
  delete [] all_energies; delete [] trial_energies;

  return 0;
}
