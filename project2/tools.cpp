#include <armadillo>
#include <math.h>
#include <string>
#include "jacobi_rotation.hpp"

void write_to_file(arma::mat &matrix, std::string filename)
{
  // Write matrix to csv file, for plotting in Python
  matrix.save(filename, arma::csv_ascii);
}

double test_tridiagonal(int n)
{
  // Test that the matrix is set up correctly, by checking eigenvalues using Arma
  double analytical_eigenval = 0;
  double mean_diff = 0;
  double a, d, h, hh;
  h = 1.0/(n*1.0);
  hh = h*h;
  a = -1.0/hh;
  d = 2.0/hh;
  arma::vec eigenvals = arma::zeros(n);
  arma::mat eigenvecs = arma::zeros(n,n);
  arma::mat tridiag = arma::zeros(n, n);
  tridiag.diag() += d;
  tridiag.diag(-1) -= a;
  tridiag.diag(+1) -= a;
  arma::eig_sym(eigenvals, eigenvecs, tridiag);
  for(int i=1; i < n+1; i++)
  {
    analytical_eigenval = d + 2*a*cos(i*M_PI/(n+1));
    mean_diff += fabs(eigenvals[i-1] - analytical_eigenval);
  }
  mean_diff = (double) mean_diff/n;
  return mean_diff;
}

double test_max_nondiagonal(int n)
{
  // Test if max_nondiagonal can find greatest non-diagonal element
  // Compare result with armadillo max function
  int max_index[2] = {-1, -1}; // save indices of max
  double max_val = 0;
  arma::mat test_matrix = arma::randn(n,n);
  arma:: mat abs_matrix = arma::abs(test_matrix);
  max_nondiagonal(test_matrix, max_index, max_val, n);
  abs_matrix.diag() = arma::zeros(n);
  double arma_result = abs_matrix.max();
  double algo_result = fabs(max_val);
  double diff = fabs(arma_result - algo_result);
  return diff;
}

double test_norm(double tolerance, int n)
{
  // Test that norm is preserved for Jacobi Rotation Algorithm.
  double norm_before, norm_after, diff, h, hh;
  tolerance = 1e-5;
  arma::mat tridiag = arma::zeros(n, n);
  h = 1.0/(n*1.0);
  hh = h*h;
  tridiag.diag() += 2.0/hh;
  tridiag.diag(-1) -= 1.0/hh;
  tridiag.diag(+1) -= 1.0/hh;
  norm_before = arma::norm(tridiag, "fro");
  jacobi_rotation(tridiag, tolerance, n);
  norm_after = arma::norm(tridiag, "fro");
  diff = fabs(norm_after - norm_before);
  return diff;
}

double test_eigenvals(double tolerance, int n)
{
  double h, hh, diff;
  arma::mat tridiag = arma::zeros(n, n);
  arma::vec eigenvals = arma::zeros(n);
  h = 1.0/(n*1.0);
  hh = h*h;
  diff = 0;
  tridiag.diag() += 2.0/hh;
  tridiag.diag(-1) -= 1.0/hh;
  tridiag.diag(+1) -= 1.0/hh;
  arma::eig_sym(eigenvals, tridiag);
  jacobi_rotation(tridiag, tolerance, n);
  tridiag.diag() = arma::sort(tridiag.diag()); // sort eigenvalues for comparison
  for(int i = 0; i < n; i++)
  {
    diff += fabs(tridiag(i,i) - eigenvals(i) );
  }
  diff = (double) diff/(n*1.0);
  return diff;
}
