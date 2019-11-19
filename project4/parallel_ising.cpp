#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include <mpi.h>
#include "ising_model.hpp"

// Parallel Ising Model, using MPI, uses Ising model class, for different temps.
// temperature range is currently hard-coded as being between 2.0 and 2.6.

// Usage: mpirun -np num_processes ./parallel_ising l n dt
// num_processes is the number of parallel processes
// l is the grid length (creates grid of quadratic lxl)
// n is the number of MC cycles/sweeps
// dt is the step in temperature

// run this for l = 40, then 60, then 80, then 100 using dt = 0.025 to plot. (takes forever)

int main(int argc, char *argv[])
{
  int num_processes = 0; // arguments for MPI, number of processes
  int process_rank = 0;  // process identifier
  int num_steps;  // number of temperature steps
  int l = atoi(argv[1]); // Size of spin grid, lxl
  int n = atoi(argv[2]); // number of Monte Carlo cycles (each of lxl spin flips)
  double dt = atof(argv[3]);
  double **all_data;
  double **final_data;
  double t; // temperature
  double t_min = 2.0;
  double t_max = 2.6;
  double time_start = 0;
  double execution_time = 0;
  // dt = (tb-ta)/(n-1) , include endpoints + 0.5f to avoid roundoff
  num_steps = (double) (t_max - t_min)/dt + 1 + 0.5f;

  // Initialize arrays/matrices to work with MPI reduce; Arma doesn't seem to work :'(
  all_data =  new double*[4]; // mean energy, mean mag, susceptibility, heat cap.
  final_data = new double*[4];  // stored in rows

  for (int i = 0; i < 4; i++) // initialize arrays in memory
  {
    all_data[i] = new double[num_steps];
    final_data[i] = new double[num_steps];
  }

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < num_steps; j++)
    {
      all_data[i][j] = 0; // initialize all elements to zero
      final_data[i][j] = 0;
    }
  }
  // init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD); // Distribute l,n, t to all procs
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  time_start = MPI_Wtime();
  // Do parallel calculation with many temperatures! Remember timing
  IsingModel model(l, n, t); // each parallel op gets its own class instance

  for(int i = 0; i < num_steps; i+= num_processes)
  {
    t = t_min + dt*(process_rank + i); // each process gets its own temperature

    if(t <= t_max ) // only calculate if temperature is in desired range
    {
      std::cout << "Process " << process_rank << " running t= " << t << std::endl;
      model.set_random_spin_config(process_rank); //different seed for each proc.
      model.set_temperature(t); // reset temperature of model
      model.metropolis(process_rank); // carry out MC simulation
      all_data[0][i + process_rank] = model.mean_energy(n-1); // save results
      all_data[1][i + process_rank] = model.mean_abs_mag(n-1);
      all_data[2][i + process_rank] = model.susceptibility(n-1);
      all_data[3][i + process_rank] = model.heat_capacity(n-1);
    }
    else
    {
      std::cout << "Process " << process_rank << " finished." << std::endl;
    }
  }

  for (int i=0; i < 4; i++)
  {   // Each process has unique elements in data matrix, sum just merges them
      MPI_Reduce(all_data[i], final_data[i], num_steps, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  execution_time = MPI_Wtime() - time_start;
  MPI_Finalize();

  if (process_rank == 0)
  {
    std::cout << "Finished in " << execution_time << " s for l=" <<l<< std::endl;
    // Write result array to file
    std::fstream filewriter;
    std::string name = "./results/parallel/"+std::to_string(l)+"x"+std::to_string(l);
    filewriter.open(name, std::ios::out); // write mode
    for(int i = 0; i < 4; i++)
    {
      for(int j = 0; j < num_steps; j++)
      {
        filewriter << final_data[i][j];
        if (j < num_steps -1)
        {
          filewriter << ","; // csv file :)
        }
      }
      filewriter << "\n";
    }
    filewriter.close();
  }

  for (int i = 0; i < 4; ++i)
  {   // clear memory
      delete [] all_data[i];
      delete [] final_data[i];
  }
  delete [] all_data; delete [] final_data;
}
