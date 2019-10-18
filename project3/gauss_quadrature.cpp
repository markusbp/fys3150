#include <iostream>
#include <cmath>
#include <armadillo>
#include <chrono>
#include <string>
#include "../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/lib.h"
#include "../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/gauss_laguerre.hpp"

void write_to_file(arma::mat &matrix, std::string filename)
{
  // Write matrix to csv file, for plotting in Python
  matrix.save(filename, arma::csv_ascii);
}

double integrand(arma::vec r1, arma::vec r2)
{
  double norm_diff = arma::norm(r1 - r2); // norm of difference
  double temp = 0;
  if (norm_diff > 1e-12) 
  {
    temp = (double) (exp(-4*(arma::norm(r1) + arma::norm(r2)))*1.0)/norm_diff;
  }
  return temp;
}

void init_zeros(double *arr, int n)
{
  for(int i = 0; i < n; i++)
  {
    arr[i] = 0;
  }
}

double cartesian_coordinates(double lower, double upper, int n)
{
  double *x, *w; // arrays for holding x-coordinates and weights, respectively
  x = new double[n];
  w = new double[n];
  init_zeros(x, n);
  init_zeros(w, n);
  gauleg(lower, upper, x, w, n);
  arma::vec r1 = arma::zeros(3);
  arma::vec r2 = arma::zeros(3);

  double integral = 0;
  for(int x1 = 0; x1 < n; x1++)
  {
    r1[0] = x[x1];
    for(int y1 = 0; y1 < n; y1++)
    {
      r1[1] = x[y1];
      for(int z1 = 0; z1 < n; z1++)
      {
        r1[2] = x[z1];
        for(int x2 = 0; x2 < n; x2++)
        {
          r2[0] = x[x2];
          for(int y2 = 0; y2 < n; y2++)
          {
            r2[1] = x[y2];
            for(int z2 = 0; z2 < n; z2++)
            {
              r2[2] = x[z2];
              integral += w[x1]*w[x2]*w[y1]*w[y2]*w[z1]*w[z2]*integrand(r1, r2);
            }
          }
        }
      }
    }
  }
  delete [] x;
  delete [] w; // clear variables
  return integral;
}

double spherical_coordinates(int n)
{
  double *x_phi, *x_theta, *x_r, *w_phi, *w_theta, *w_r; // arrays for holding x-coordinates and weights, respectively
  w_r = new double[n+1]; x_r = new double[n+1];
  w_phi = new double[n]; x_phi = new double[n];
  w_theta = new double[n]; x_theta = new double[n];

  init_zeros(w_r, n+1);
  init_zeros(x_r, n+1);
  init_zeros(w_phi, n);
  init_zeros(x_phi, n);
  init_zeros(w_theta, n);
  init_zeros(x_theta, n);

  double r1, r2, r_12, beta, gamma, theta1, theta2, phi1, phi2, sinsin;
  r1 = 0; r2 = 0; r_12 = 0; theta1 = 0; theta2 = 0;
  phi1 = 0; phi2 = 0; sinsin = 0; beta = 0, gamma = 0;

  gauss_laguerre(x_r, w_r, n, 0);
  gauleg(0, 2*M_PI, x_phi, w_phi, n);
  gauleg(0, M_PI, x_theta, w_theta, n);

  double integral = 0;
  for(int i = 0; i < n; i++)
  {
    phi1 = x_phi[i];
    for(int j = 0; j < n; j++)
    {
      phi2 = x_phi[j];
      for(int k = 0; k < n; k++)
      {
        theta1 = x_theta[k];
        for(int l = 0; l < n; l++)
        {
          theta2 = x_theta[l];
          for(int m = 0; m < n; m++)
          {
            r1 = x_r[m+1];
            for(int p = 0; p < n; p++)
            {
              r2 = x_r[p+1];
              sinsin = sin(theta1)*sin(theta2);

              beta = cos(theta1)*cos(theta2) + sinsin*cos(phi1 - phi2);
              gamma = r1*r1 + r2*r2 - 2*r1*r2*beta;
              if(gamma > 1e-12)
              {
                r_12 = sqrt(gamma);
                integral  += (double) 1.0/r_12 *exp(-3*(r1 + r2))*r1*r1*r2*r2*sinsin*w_r[p+1]*w_r[m+1]*w_theta[l]*w_theta[k]*w_phi[j]*w_phi[i];
              }
            }
          }
        }
      }
    }
  }

  delete [] x_phi; delete [] x_theta; delete [] x_r;
  delete [] w_phi; delete [] w_theta; delete [] w_r;
  return integral;
}

int main(int argc, char *argv[])
{
  int n; // number of integration points
  double lower, upper, integral; // limits of integration
  std::string integral_type;

  integral = 0;

  integral_type = argv[1];
  n = atoi(argv[2]);

  if (integral_type == "cartesian")
  {
    lower = atof(argv[3]);
    upper = atof(argv[4]);
    integral = cartesian_coordinates(lower, upper, n);
  }

  if (integral_type == "spherical")
  {
    integral = spherical_coordinates(n);
  }

  if (integral_type == "profile_legendre" || integral_type == "profile_laguerre")
  {
    int start = 2;
    int num_steps = n;
    int step = 4;
    lower = atof(argv[3]);
    upper = atof(argv[4]);
    arma::vec results = arma::zeros(num_steps); // integration results for different n
    arma::vec params = arma::zeros(num_steps + 2); // Parameters; limits, start and stop
    params[0] = lower; params[1] = upper;

    if (integral_type == "profile_legendre")
    {
      for(int i = 0; i < num_steps; i++)
      {
        int index = i*step + start;
        results[i] = cartesian_coordinates(lower, upper, index);
        params[i+2] = index;
        std::cout << index << std::endl;
      }
      write_to_file(results, "legendre_profile_results");
      write_to_file(params, "legendre_profile_params");
    }
    else
    {
      for(int i = 0; i < num_steps; i++)
      {
        int index = i*step + start;
        results[i] = spherical_coordinates(index);
        params[i+2] = index;
        std::cout << index << std::endl;
      }
      write_to_file(results, "laguerre_profile_results");
      write_to_file(params, "laguerre_profile_params");
    }
  }


  if (integral_type == "spherical" || integral_type == "cartesian")
  {
    std::cout << "Integral is approximately : " << integral << std::endl;
  }
  return 0;
}
