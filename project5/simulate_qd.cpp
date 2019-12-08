#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include <mpi.h>
#include "quantum_dots.hpp"

int main(int argc, char *argv[])
{
  int process_rank = 0;
  int num_processes = 0;
  double time_start = 0;
  double execution_time = 0;

  int n = 5000000;
  int col = 0;
  int grid_length = 100;
  double freq = 1;
  double grid_step = 0.02;
  double beta_init = 0.01;
  double alpha_init = 0.01;
  arma::mat variational_params = arma::zeros(2, grid_length);

  // Initialize arrays/matrices to work with MPI reduce
  double **all_energies;
  double **trial_energies;
  all_energies =  new double*[grid_length];
  trial_energies =  new double*[grid_length];

  for (int i = 0; i < grid_length; i++)
  {
    all_energies[i] = new double[grid_length];
    trial_energies[i] = new double[grid_length];
    variational_params(0, i) = (i+1)*grid_step;
    variational_params(1, i) = (i+1)*grid_step;
    for (int j = 0; j < grid_length; j++)
    {
      trial_energies[i][j] = 0; // initialize all elements to zero
      all_energies[i][j] = 0; // initialize all elements to zero
    }
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  time_start = MPI_Wtime();

  Psi2 trial2(n, alpha_init, beta_init, freq, process_rank);

  for(int i = 0; i < grid_length; i++)
  {
    trial2.a = variational_params(0, i);
    trial2.varparam = beta_init;
    for(int j = 0; j < grid_length; j+= num_processes)
    {
      col = j + process_rank;
      if(i < grid_length && col < grid_length)
      {
        trial2.varparam = variational_params(1, col);
        trial2.metropolis();
        trial_energies[i][col] = trial2.averages(n, 0);
      }
    }
    if(process_rank==0)
    {
      std::cout << "Computing row " << i+1 << "/" << grid_length << std::endl;
    }
  }

  for (int i = 0; i < grid_length; i++)
  { // Each process has unique elements in data matrix, sum just merges them
    MPI_Reduce(trial_energies[i], all_energies[i], grid_length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  execution_time = MPI_Wtime() - time_start;

  MPI_Finalize();

  if(process_rank == 0)
  {
    // Write result array to file
    std::fstream filewriter;
    std::string name = "./results/full_gridsearch_psi2";
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
    variational_params.save("results/full_gs_psi2_varparams", arma::csv_ascii);
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
