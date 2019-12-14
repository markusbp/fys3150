#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]); // mc cycles
  int mc_cycles = atoi(argv[2]);
  double alpha = 1.05; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double beta = 0.23; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double start_freq = 0.01;
  double stop_freq = 1.00;
  double freq_step = (double)(stop_freq - start_freq)/n;
  Psi2 trial2(mc_cycles, alpha, beta, start_freq, 0); // n, alpha, omega, seed
  arma::mat params(4, n); // potential energy, kin. energy, values of omega

  for(int j = 0; j < n; j++)
  {
    params(3, j) = (j+1)*freq_step; // Different values of Omega
  }

  for(int i = 0; i < n; i++)
  {
    trial2.set_omega(params(3,i));
    trial2.metropolis();
    params(0, i) = trial2.averages(mc_cycles, 2); // Ek
    params(1, i) = trial2.averages(mc_cycles, 3); // Ep, interacting
    params(2, i) = trial2.averages(mc_cycles, 6); // Ep, non-interacting
  }

  params.save("results/psi2/virial", arma::csv_ascii);
  return 0;
}
