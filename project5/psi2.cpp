#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]); // mc cycles
  double freq = atof(argv[2]); // harmonic oscillator frequency
  double alpha = atof(argv[3]); // initial value
  double beta = atof(argv[4]);
  // Instantiate Quantum Dot model with Psi1 trial wave function
  Psi2 trial2(n, alpha, beta, freq, 0); // n, alpha, omega, seed
  trial2.metropolis();

  trial2.averages.save("results/psi2_run", arma::csv_ascii);
  return 0;
}
