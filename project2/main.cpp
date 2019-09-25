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
  num_vals = 20;
  arma::vec duration = arma:: zeros(num_vals);
  for(int i = 0; i < num_vals; i++)
  {
    h = (double) 1.0/(n+1);
    hh = h*h;
    arma::mat tridiag = arma::zeros(n, n);
    tridiag.diag() += 2.0/hh;
    tridiag.diag(-1) -= 1.0/hh;
    tridiag.diag(+1) -= 1.0/hh;
    auto start = std::chrono::high_resolution_clock::now(); // time execution
    jacobi_rotation(tridiag, tolerance, n);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    duration(i) = elapsed.count();
    tridiag.zeros();
    n = n+10;
  }
  std::cout << duration << std::endl;
  write_to_file(duration, "timedata");
}

int main(int argc, char *argv[])
{
  int n;
  double h, hh, r_max, r0;
  r_max = 10.0;
  r0 = 0.0;
  n = atoi(argv[1]);
  h = (double) (r_max - r0)/n;
  hh = h*h;
  n = (int) n-1; // number of points for which calculation is actually performed
  // Initialize all elements of rho (as vector).
  arma::vec r = arma::zeros(n);
  for (int i = 0; i < n ; i++)
  {
    r[i] = (i+1.0)*h;
  }

  arma::vec eigenvals = arma::zeros(n);
  arma::mat eigenvecs = arma::zeros(n,n);
  arma::mat tridiag = arma::zeros(n, n);

  // With potential rho^2
  tridiag.diag() += 2.0/(hh) + arma::square(r);
  tridiag.diag(-1) -= 1.0/hh;
  tridiag.diag(+1) -= 1.0/hh;

  double a = -1.0/hh;
  double d = 2.0/hh;

  std::cout << "Norm before: " << norm(tridiag, "fro") << std::endl;
  double eps = 1.0e-5; // tolerance
  jacobi_rotation(tridiag, eps, n);
  std::cout<< "Norm after: " << norm(tridiag, "fro") << std::endl;
  std::cout << "Diagonal elements (sorted) \n" << arma::sort(tridiag.diag()) << std::endl;
  write_to_file(tridiag, "test");

  time_jacobi_rotation();
  return 0;

}
