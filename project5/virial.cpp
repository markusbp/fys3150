#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"
// Script to compute and save kinetic and potential energy for optimal psi2 wave func.
// Used to generate data for testing virial theorem for different HO frequencies

// usage:
// ./virial num_omegas mc_cycles
// num_omegas: number of different omega values between 0.01 and 1.00 to be used
// mc_cycles: Number of Monte Carlo simulation cycles. Usually 10^6.
int main(int argc, char *argv[])
{
  int n = atoi(argv[1]);
  int mc_cycles = atoi(argv[2]);

  double alpha1 = 1.00;

  double alpha2 = 1.000; // (empirical) optimal alpha and beta
  double beta2 = 0.270;

  double start_freq = 0.01; // starting omega
  double stop_freq = 1.00;  // maximum omega
  double freq_step = (double)(stop_freq - start_freq)/n;
  // Instantiate Psi2 wavefunction quantum dot model
  Psi2 trial2(mc_cycles, alpha2, beta2, start_freq, 0); // n, alpha, omega, seed
  arma::mat params(7, n); // (kin. energy, potential energy, non.int.)*2, values of omega

  // Instantiate Psi1 wavefunction quantum dot model
  Psi1 trial1(mc_cycles, alpha1, start_freq, 0); // n, alpha, omega, seed

  for(int j = 0; j < n; j++)
  {
    params(3, j) = (j+1)*freq_step; // Different values of Omega
  }

  for(int i = 0; i < n; i++)
  {
    trial2.set_omega(params(3,i)); // setter method for changing omega
    trial2.metropolis();
    params(0, i) = trial2.averages(mc_cycles, 2); // Ek
    params(1, i) = trial2.averages(mc_cycles, 3); // Ep, interacting
    params(2, i) = trial2.averages(mc_cycles, 6); // Ep, non-interacting

    trial1.set_omega(params(3,i)); // setter method for changing omega
    trial1.metropolis();
    params(4, i) = trial1.averages(mc_cycles, 2); // Ek
    params(5, i) = trial1.averages(mc_cycles, 3); // Ep, interacting
    params(6, i) = trial1.averages(mc_cycles, 6); // Ep, non-interacting
  }
  params.save("results/psi2/virial", arma::csv_ascii);
  return 0;
}
