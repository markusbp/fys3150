#include <iostream>
#include <iomanip>
#include <armadillo>
#include <stdlib.h>
#include <math.h>
#include "tools.hpp"
#include "jacobi_rotation.hpp"

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
  write_to_file(tridiag, "test", n);
  return 0;

}
