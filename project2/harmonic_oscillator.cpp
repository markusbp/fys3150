#include <iostream>
#include <iomanip>
#include <armadillo>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <string>
#include "tools.hpp"
#include "jacobi_rotation.hpp"
  
/* Multi-purpose main function for calculating eigenvalues and normalized
eigenstates of one- and two-electron quantum energy eigenvalue problems.
Also times Jacobi's method implementation, and finds the mean error for
eigenvalue calculations using Jacobi's method for a grid of rho_max, n values.

Usage:

./harmonic_oscillator n  rho_max program_type -omega_r

- n is number of calculation steps
- rho_max is approximation of infinity
- program_type is the type of program to run, can be either "1particle",
  "2particle", "timing", or "error_surface".
- omega_r is optional argument, only needed for "2particle" case, and reflects
  strength of Coloumbic interaction between particles.

Example run:
./harmonic_oscillator 350 8 2particle 0.01
-- runs two-particle simulation, with n = 350 and rho_max = 8, with omega_r = 0.01
Saves result of eigenvectors, rho and eigenvalues as csv files in current directory.
*/

void time_jacobi_rotation()
{
  // Time the jacobi rotation algorithm for different matrix sizes
  double h, hh, tolerance;
  tolerance = 1e-10; // Set tolerance
  int n, num_vals;
  n = 5; // initial n value
  num_vals = 20; // number of values to test
  // initialize matrices and vectors
  arma::vec duration = arma:: zeros(num_vals);
  arma::vec n_values = arma:: zeros(num_vals); // for saving
  arma::mat tridiag = arma::zeros(n, n);
  std::cout << "Starting timing..." << std::endl;
  for(int i = 0; i < num_vals; i++)
  {
    std::cout << "n = " << n << std::endl;
    h = (double) 1.0/(n+1);
    hh = h*h;
    tridiag.diag() += 2.0/hh;
    tridiag.diag(-1) -= 1.0/hh;
    tridiag.diag(+1) -= 1.0/hh; // Initialize non-diagonal elements
    auto start = std::chrono::high_resolution_clock::now(); // time execution
    jacobi_rotation(tridiag, tolerance, n);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    duration(i) = elapsed.count(); // save execution time for a given n
    n_values(i) = n;
    n = n+20; // advance n
    tridiag.zeros(n,n); // set to zero AND resize tridiagonal
  }
  write_to_file(duration, "timedata"); // Save results and parameters
  write_to_file(n_values, "n_values_timing");
}

void find_error_surface()
{
  // Find error for different values of matrix size n and infinity approx. rhomax
  double h, hh, tolerance, rho_max, initial_rho, rho_step;
  tolerance = 1e-10; // Set tolerance
  int n, n_step, num_vals, eigenvals_to_mean;

  num_vals = 11;  // check num_valsxnum_vals grid of values

  // Initialize variables
  arma::mat rho_max_values = arma::zeros(num_vals, num_vals);
  arma::mat mean_error = arma::zeros(num_vals, num_vals);
  arma::mat n_values = arma::zeros(num_vals, num_vals);

  n_step = 35; // stepsize, n
  n = 5; // initial n value
  arma::mat tridiag = arma::zeros(n, n);

  rho_step = 2; // stepsize, rho_max
  initial_rho = 1; // initial rho value
  rho_max = initial_rho;
  // Find mean error for the first eigenvals_to_mean eigenvalues.
  eigenvals_to_mean = 1; // must be less than or equal to initial n value
  arma::vec rho_vec = arma::zeros(n);
  std::cout << "Started finding error surface..." << std::endl;

  for(int i = 0; i < num_vals; i++)
  {
    std::cout << "n = "<< n << std::endl;
    for(int j = 0; j < num_vals; j++)
    {
      tridiag.zeros(n,n); // Set all elements to zero and resize
      rho_vec.zeros(n);
      h = (double) rho_max/(n*1.0+1.0);   // rho_0 zero
      hh = h*h;
      // Initialize rho vector for new n
      for (int l = 0; l < n ; l++)
      {
        rho_vec[l] = (l+1.0)*h;
      }
      tridiag.diag() += 2.0/(hh*1.0) + arma::square(rho_vec);
      tridiag.diag(-1) -= 1.0/(hh*1.0);
      tridiag.diag(+1) -= 1.0/(hh*1.0);
      jacobi_rotation(tridiag, tolerance, n); // diagonalize
      tridiag.diag() = arma::sort(tridiag.diag()); // sort in ascending order
      n_values(i, j) = n; // save n value
      rho_max_values(i, j) = rho_max; // save rho value
      for(int k = 0; k < eigenvals_to_mean; k++)
      {
        // find mean (absolute) error
        mean_error(i, j) += fabs(tridiag(k, k) - (4*k + 3)); // Analytical eigenvalue
      }
      mean_error(i, j) = (double) mean_error(i, j)/(eigenvals_to_mean*1.0);
      rho_max += rho_step;
    }
    n += n_step;
    rho_max = initial_rho;
  }
  write_to_file(rho_max_values, "max_rho_values"); // save parameters and error
  write_to_file(n_values, "n_values_error_surface");
  write_to_file(mean_error, "mean_eigenvalue_error");
}

void normalize_eigenstates(arma::mat &eigenstate_matrix, arma::vec &rho, int n)
{
  // Find normalized eigenstates over interval defined by rho
  // Compute norm from integral
  arma::mat norm = arma::trapz(rho, arma::square(eigenstate_matrix), 0);
  for(int i=0; i < n; i++)
  {
    // eigenstates are saved in columns, normalize by dividing by norm
    eigenstate_matrix.col(i) = eigenstate_matrix.col(i)*(1.0/sqrt(norm(0, i)));
  }
}

int main(int argc, char *argv[])
{
  int n;
  double h, hh, r_max, r0, tolerance;
  std::string program_type;

  tolerance = 1e-10; // Error tolerance for jacobi rotation
  r0 = 0.0; // origin
  n = atoi(argv[1]); // Number of calculation points
  r_max = atof(argv[2]); // Approximation of infinity
  program_type = argv[3]; // Can be 1d, 2d, timing, or error_surface

  h = (double) (r_max - r0)/(n*1.0 + 1.0); // stepsize
  hh = h*h;
  // Initialize all elements of rho (as vector), for saving.
  arma::vec rho = arma::zeros(n);
  for (int i = 0; i < n ; i++)
  {
    rho[i] = (i+1.0)*h;
  }
  // Initialize tridiagonal
  arma::vec eigenvals = arma::zeros(n);
  arma::mat eigenvecs = arma::zeros(n,n);
  arma::mat tridiag = arma::zeros(n, n);

  if (program_type == "1particle")
  {
    // 3D quantum dot, one electron case
    tridiag.diag() += 2.0/(hh) + arma::square(rho); // Potential
    tridiag.diag(-1) -= 1.0/hh;
    tridiag.diag(+1) -= 1.0/hh;
    arma::eig_sym(eigenvals, eigenvecs, tridiag);
    normalize_eigenstates(eigenvecs, rho, n);
    jacobi_rotation(tridiag, tolerance, n);
    arma::mat diagonal_elements = arma::sort(tridiag.diag());
    write_to_file(rho, "eigenvecs_rho_1");
    write_to_file(eigenvecs, "eigenvec_1");
    write_to_file(diagonal_elements, "jacobi_eigenvals_1");
  }

  else if (program_type == "2particle")
  {
    double omega = atof(argv[4]);
    // 3D quantum dot, two electron case
    // Potential
    tridiag.diag() += 2.0/hh + omega*omega*arma::square(rho) + arma::ones(n)/rho;
    tridiag.diag(-1) -= 1.0/hh;
    tridiag.diag(+1) -= 1.0/hh;
    // Arma solution
    arma::eig_sym(eigenvals, eigenvecs, tridiag);
    normalize_eigenstates(eigenvecs, rho, n);
    // Jacobi's method for finding eigenvalues
    jacobi_rotation(tridiag, tolerance, n);
    arma::mat diagonal_elements = arma::sort(tridiag.diag());
    // Give name useful for plotting
    std::string eigenvec_name = "eigenvecs_2_omega_";
    std::string eigenval_name = "eigenvals_2_omega_";
    eigenvec_name.append(argv[4]);
    eigenval_name.append(argv[4]);
    write_to_file(rho, "eigenvecs_rho_2");
    write_to_file(eigenvecs, eigenvec_name);
    write_to_file(diagonal_elements, eigenval_name);
  }
  else if (program_type == "timing"){time_jacobi_rotation();} // time execution
  else if (program_type == "error_surface"){find_error_surface();}
  return 0;
}
