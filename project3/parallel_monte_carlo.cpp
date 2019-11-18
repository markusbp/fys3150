#include <cmath>
#include <mpi.h>
#include <iostream>
#include <armadillo>

/*
  Script for running multiple Monte Carlo simulations in parallel,
  calculating expected value of particle separation in Helium atom, an integral.
  Usage:
  mpiexec -np num_processes ./parallel_monte_carlo n
  - num_processes is the number of parallel processes, I use 8
  - n is the number of integration points, in total. Must be divisible by num_processes
  The script appends the value for the integral to a file. To be used in
  conjuction with profile_parallel_mc bash script for profiling!
*/

int main(int argc, char *argv[])
{
  std::string program_type;

  int n = atoi(argv[1]);
  double execution_time = 0;
  double integral = 0;

  int num_processes = 0;
  int process_rank  = 0;
  double start = 0, stop = 0;
  double r1 = 0, r2 = 0, r1r1 = 0, r2r2 = 0, r_12 = 0, alpha = 0;
  double theta1 = 0, theta2 = 0, phi1 = 0, phi2 = 0, sinsin = 0;
  double local_integral = 0, total_integral = 0;

  double jacobi_determinant = 4*pow(M_PI, 4);
  double TWO_PI = 2*M_PI;

  // init MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // All processes share n
  srand(time(NULL) + process_rank); // Initialize RNG
  int n_per_process = n/num_processes;

  if(n%num_processes!=0 && process_rank == 0)
  {
    // demand that processes are distributed equally
    std::cout << "Number of processes and integration points not compatible," <<
                  "must be equal number of calculations per process." << std::endl;
    exit(1);
  }
  start = MPI_Wtime(); // capture execution time
  for(int i = 0; i < n_per_process; i++)
  {
    // Same calculation as for importance algorithm
    r1 = - log( 1 - (double) rand()/RAND_MAX); // create random numbers, avoid user-defined funcs.
    r2 = - log( 1 - (double) rand()/RAND_MAX);
    theta1 = (double) M_PI*rand()/(RAND_MAX);
    theta2 = (double) M_PI*rand()/(RAND_MAX);
    phi1 = (double) TWO_PI*rand()/(RAND_MAX);
    phi2 = (double) TWO_PI*rand()/(RAND_MAX);
    r1r1 = r1*r1; // precalc to save computation
    r2r2 = r2*r2;
    sinsin = sin(theta1)*sin(theta2);
    alpha  = cos(theta1)*cos(theta2) + sinsin*cos(phi1 - phi2);
    r_12 = sqrt(r1r1 + r2r2 - 2*r1*r2*alpha);
    local_integral  += (double) 1.0/r_12*r1r1*r2r2*sinsin*exp(-3*(r1 + r2));
  }
  // merge all local integrals to total integral by summing each contribution
  MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total_integral = (double) total_integral/(n)*jacobi_determinant;
  stop = MPI_Wtime();

  if(process_rank == 0) // only print if master node
  {
    // Append result to file, works in conjunction with bash script profile_parallel_mc
    integral = total_integral;
    execution_time = stop - start;
    std::cout << "n = " <<  n << std::endl;
    std::cout << "Total integral is approximately " << integral << std::endl;
    std::cout << "Computation finished in " << execution_time << std::endl;

    std::ofstream save_integral, save_timing, save_n;
    save_integral.open ("results/parallel_monte_carlo_integral", std::ios::app);
    save_timing.open ("results/parallel_monte_carlo_timing", std::ios::app);
    save_n.open ("results/parallel_monte_carlo_n_values", std::ios::app);
    save_integral << integral;
    save_timing << execution_time;
    save_n << n;
    save_timing.close(); save_integral.close(), save_n.close();
  }
  MPI_Finalize(); // end simulation

  return 0;
}
