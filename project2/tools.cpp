#include <armadillo>
#include <math.h>
#include <string>
#include "jacobi_rotation.hpp"

void write_to_file(arma::mat &matrix, std::string filename, int n)
{
  // Write matrix to binary file, for plotting in Python
  matrix.save(filename, arma::raw_binary);
}

double test_eigenval_decomp(arma::mat &sym_matrix, arma::vec &eigenvals,
                          arma::mat &eigenvecs, double a, double d, int n)
{
  // Test that the matrix is set up correctly, by checking eigenvalues using Arma
  double analytical_eigenval = 0;
  double mean_diff = 0;
  arma::eig_sym(eigenvals, eigenvecs, sym_matrix);
  for(int i=1; i < n+1; i++)
  {
    analytical_eigenval = d + 2*a*cos(i*M_PI/(n+1));
    mean_diff += fabs(eigenvals[i-1] - analytical_eigenval);
  }
  mean_diff = (double) mean_diff/n;
  return mean_diff;
}

void test_max_nondiagonal(int n)
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
  std::cout << diff << std::endl;
  std::cout << max_val << std::endl;
}

double test_norm(double tolerance, int n)
{
  // Test that norm is preserved for Jacobi Rotation Algorithm.
  double norm_before, norm_after, diff;
  tolerance = 1e-5;
  arma::mat tridiag = arma::zeros(n, n);
  tridiag.diag() += 2.0;
  tridiag.diag(-1) -= 1.0;
  tridiag.diag(+1) -= 1.0;
  norm_before = arma::norm(tridiag, "fro");
  jacobi_rotation(tridiag, tolerance, n);
  norm_after = arma::norm(tridiag, "fro");
  diff = fabs(norm_after - norm_before);
  return diff;
}
