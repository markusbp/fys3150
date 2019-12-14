#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"

double grid_search_psi1(Psi1 instance, double start, double nudge, int steps, std::string name)
{
  double min_en  = 1e6; // Random large initial value
  double curr_en = 0; // current energy
  double min_par = 0; // parameter corresponding to minimal energy
  double n = instance.mc_cycles; // number of MC cycles
  instance.reset(); // clear all matrices just to be safe :)
  instance.varparam = start; // set initial value of alpha

  arma::mat grid_params = arma::zeros(3, steps); // save energy and alpha

  for(int i = 0; i<steps; i++)
  {
    // run MC simulation with metropolis algorithm
    instance.metropolis();
    curr_en = instance.averages(n, 0); // get energy
    if(curr_en < min_en)
    {
      // if smaller energy, just for printing
      min_en = curr_en;
      min_par = instance.varparam;
    }
    grid_params(0, i) = curr_en;
    grid_params(1, i) = instance.varparam; // save quantities
    grid_params(2, i) = instance.averages(n, 1);
    instance.reset(); // clear for new iteration
    // reset average and positions
    instance.varparam += nudge; // change variational parameter slightly
  }
  std::cout << "Alpha "<<  min_par <<" Min. Energy "<< min_en << std::endl;
  grid_params.save(name, arma::csv_ascii); // save all
  return min_par;
}

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]); // mc cycles
  double freq = atof(argv[2]); // harmonic oscillator frequency
  double alpha_init = atof(argv[3]); // initial value
  double best_alpha = 0;
  // Instantiate Quantum Dot model with Psi1 trial wave function
  Psi1 trial1(n, alpha_init, freq, 0); // n, alpha, omega, seed

  std::string freqstring(argv[2]);
  std::string name = "results/psi1/psi1_gridsearch_w=" + freqstring;
  best_alpha = grid_search_psi1(trial1, alpha_init, 0.01, 500, name); // perform search for best alpha

  Psi1 run2(n, best_alpha, freq, 0);
  run2.metropolis();
  run2.averages.save("results/psi1/psi1_opt_w=" + freqstring, arma::csv_ascii);
  return 0;
}
