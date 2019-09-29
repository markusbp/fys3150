#include <iostream>
#include <iomanip>
#include <armadillo>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include "tools.hpp"
#include "jacobi_rotation.hpp"

void time_jacobi_rotation()
{
  double h, hh, tolerance;
  tolerance = 1e-5;
  int n, num_vals;
  n = 5;
  num_vals = 25;
  arma::vec duration = arma:: zeros(num_vals);
  arma::vec n_values = arma:: zeros(num_vals);
  arma::mat tridiag = arma::zeros(n, n);
  std::cout << "Starting timing..." << std::endl;
  for(int i = 0; i < num_vals; i++)
  {
    std::cout << "n = " << n << std::endl;
    h = (double) 1.0/(n+1);
    hh = h*h;
    tridiag.diag() += 2.0/hh;
    tridiag.diag(-1) -= 1.0/hh;
    tridiag.diag(+1) -= 1.0/hh;
    auto start = std::chrono::high_resolution_clock::now(); // time execution
    jacobi_rotation(tridiag, tolerance, n);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    duration(i) = elapsed.count();
    n_values(i) = n;
    n = n+10;
    tridiag.zeros(n,n);
  }
  write_to_file(duration, "timedata");
  write_to_file(n_values, "n_values_timing");
}

void find_error_surface()
{
  double h, hh, tolerance, rho_max, initial_rho, rho_step;
  tolerance = 1e-5;
  int n, n_step, num_vals, eigenvals_to_mean;

  num_vals = 15;

  arma::mat rho_max_values = arma::zeros(num_vals, num_vals);
  arma::mat mean_error = arma::zeros(num_vals, num_vals);
  arma::mat n_values = arma::zeros(num_vals, num_vals);

  n_step = 10;
  n = 5; // initial n value
  arma::mat tridiag = arma::zeros(n, n);

  rho_step = 1;
  initial_rho = 1;
  rho_max = initial_rho;

  eigenvals_to_mean = 1; // must be less than or equal to initial n value
  arma::vec analytical_eigenvals = {3};
  std::cout << "Started finding error surface..." << std::endl;

  for(int i = 0; i < num_vals; i++)
  {
    for(int j = 0; j < num_vals; j++)
    {
      tridiag.zeros(n,n); // Set all elements to zero
      h = (double) rho_max/(n+1);   // rho_0 zero
      hh = h*h;
      tridiag.diag() += 2.0/hh;
      tridiag.diag(-1) -= 1.0/hh;
      tridiag.diag(+1) -= 1.0/hh;
      jacobi_rotation(tridiag, tolerance, n);
      n_values(i, j) = n;
      rho_max_values(i, j) = rho_max;

      for(int k = 0; k < eigenvals_to_mean; k++)
      {
        mean_error(i, j) += fabs(tridiag(k, k) - analytical_eigenvals(k));
      }
      mean_error(i, j) = (double) mean_error(i, j)/eigenvals_to_mean;

      rho_max += rho_step;
    }
    n += n_step;
    tridiag.zeros(n, n); // Set all elements to zero AND resize
    rho_max = initial_rho;
    std::cout << "n = "<< n << std::endl;
  }
  write_to_file(rho_max_values, "max_rho_values");
  write_to_file(n_values, "n_values_error_surface");
  write_to_file(mean_error, "mean_eigenvalue_error");
}

void normalize_eigenstates(arma::mat &eigenstate_matrix, arma::vec &rho, int n)
{
  std::cout<< size(eigenstate_matrix)<< size(rho) << std::endl;
  arma::mat norm = arma::trapz(rho, arma::square(eigenstate_matrix), 0);
  for(int i=0; i < n; i++)
  {
    eigenstate_matrix.col(i) = eigenstate_matrix.col(i)*(1.0/sqrt(norm(0, i)));
  }

}

int main(int argc, char *argv[])
{
  int n;
  double h, hh, r_max, r0;
  r_max = 8.0;
  r0 = 0.0;
  n = atoi(argv[1]);
  h = (double) (r_max - r0)/n;
  hh = h*h;
  n = (int) n-1; // number of points for which calculation is actually performed
  // Initialize all elements of rho (as vector).
  arma::vec rho = arma::zeros(n);
  for (int i = 0; i < n ; i++)
  {
    rho[i] = (i+1.0)*h;
  }
  // Initialize tridiagonal
  arma::vec eigenvals = arma::zeros(n);
  arma::mat eigenvecs = arma::zeros(n,n);
  arma::mat tridiag = arma::zeros(n, n);

  /*
  // 3D quantum dot, one electron case
  tridiag.diag() += 2.0/(hh) + arma::square(rho);
  tridiag.diag(-1) -= 1.0/hh;
  tridiag.diag(+1) -= 1.0/hh;
  */

  double omega = 0.5;
  // 3D quantum dot, two electron case
  tridiag.diag() += omega*omega*arma::square(rho) + arma::ones(n)/rho;
  tridiag.diag(-1) -= 1.0/hh;
  tridiag.diag(+1) -= 1.0/hh;

  arma::eig_sym(eigenvals, eigenvecs, tridiag);
  std::cout<< eigenvecs << std::endl;
  normalize_eigenstates(eigenvecs, rho, n);
  std::cout<< eigenvecs << std::endl;
  //std::cout<< eigenvals << std::endl;
  // Misc.
  //time_jacobi_rotation();
  //find_error_surface();
  write_to_file(eigenvecs, "eigenvectors_2d");
  write_to_file(rho, "eigenvecs_rho");
  return 0;
}
