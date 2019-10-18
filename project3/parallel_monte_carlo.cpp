#include <cmath>
#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[])
{

  int num_processes = 0;
  int process_rank  = 0;
  double r1 = 0, r2 = 0, r1r1 = 0, r2r2 = 0, r_12 = 0, alpha = 0;
  double theta1 = 0, theta2 = 0, phi1 = 0, phi2 = 0, sinsin = 0;
  double local_integral = 0, total_integral = 0;

  double jacobi_determinant = 4*pow(M_PI, 4);
  double TWO_PI = 2*M_PI;

  int n = atoi(argv[1]);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  srand(time(NULL) + process_rank); // Initialize RNG

  for(int i = 0; i < n; i++)
  {
    r1 = - log( 1 - (double) rand()/RAND_MAX);
    r2 = - log( 1 - (double) rand()/RAND_MAX);

    theta1 = (double) M_PI*rand()/(RAND_MAX);
    theta2 = (double) M_PI*rand()/(RAND_MAX);
    phi1 = (double) TWO_PI*rand()/(RAND_MAX);
    phi2 = (double) TWO_PI*rand()/(RAND_MAX);

    r1r1 = r1*r1;
    r2r2 = r2*r2;
    sinsin = sin(theta1)*sin(theta2);
    alpha  = cos(theta1)*cos(theta2) + sinsin*cos(phi1 - phi2);
    r_12 = sqrt(r1r1 + r2r2 - 2*r1*r2*alpha);
    local_integral  += (double) 1.0/r_12*r1r1*r2r2*sinsin*exp(-3*(r1 + r2));
  }

  MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total_integral = (double) total_integral/(n*num_processes)*jacobi_determinant;
  if(process_rank == 0)
  {
    std::cout << "Total integral is approximately " << total_integral << std::endl;
  }
  MPI_Finalize();
  return 0;

}
