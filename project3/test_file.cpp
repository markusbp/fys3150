#include <cmath>
#include "../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/lib.h"
#include "../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/gauss_laguerre.hpp"

// Test functions, for unit tests. Mainly checks understanding, and usage of library functions.

double test_gauss_legendre(int n)
{
  // Test understanding of gauss_legendre/sanity of library function

  double *x, *w; // arrays for holding x-coordinates and weights, respectively
  x = new double[n];
  w = new double[n];
  double integral = 0;
  double exact_val = 0;

  double lower = -1;
  double upper = 1;
  exact_val = 2; // integral from -1 to 1 of (x+1) = 2

  for(int i = 0; i < n; i++)
  {
    x[i] = 0;
    w[i] = 0;
  }
  gauleg(lower, upper, x, w, n);

  for(int i = 0; i < n; i++)
  {
    integral += w[i]*(x[i] + 1); // function x + 1
  }

  return fabs(integral - exact_val);
}

double test_gauss_laguerre(int n)
{
  // Test understanding of gauss_laguerre/sanity of library function
  double *x, *w; // arrays for holding x-coordinates and weights, respectively
  x = new double[n+1];
  w = new double[n+1];
  double integral = 0;
  double exact_val = 0;

  exact_val = 1 ; // integral from 0 to infinity of e^(-x) = 1

  for(int i = 0; i <= n; i++)
  {
    x[i] = 0;
    w[i] = 0;
  }
  gauss_laguerre(x, w, n, 0);

  for(int i = 0; i < n; i++)
  {
    integral += w[i+1]; // function x + 1
  }
  return fabs(exact_val - integral);
}

double test_monte_carlo(int n)
{
  // test simple monte carlo integration to verify implementation
  double integral  = 0;
  double exact_val = 1.0;
  double x = 0;
  srand(time(NULL)); // Initialize RNG

  for(int i = 0; i < n; i++)
  {
    x = (double) rand()/RAND_MAX;
    integral += 4*pow(x, 3);
  }
  integral = (double) integral/n;

  return fabs(integral - exact_val);

}
