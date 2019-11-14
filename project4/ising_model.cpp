#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "ising_model.hpp"

IsingModel::IsingModel(int l, int n, double t):
            spins(l, l), // initialize fixed-size armadillo vectors/matrices
            transition_probabilites(17),
            susceptibility(n+1),
            heat_capacity(n+1),
            mean_mag(n+1),
            mean_abs_mag(n+1),
            mean_energy(n+1),
            all_energies(n+1),
            accepted_flips(n+1)
            // Constructor
{
  // Initialize relevant parameters
  kb =  1.0; // Boltzmann Constant, scaled
  temperature = t;
  beta = (double) 1.0/(kb*temperature);
  mc_cycles = n;
  lattice_size = l;
  gridsize = lattice_size*lattice_size;
  set_ordered_spin_config("all_up"); // Default configuration
}

void IsingModel::initialize_parameters()
{
  susceptibility.fill(0);
  heat_capacity.fill(0);
  mean_mag.fill(0);
  mean_abs_mag.fill(0);
  mean_energy.fill(0);
  accepted_flips.fill(0);
  energy = 0;
  magnetization = 0;
  flips = 0;

  for(int y = 0; y < lattice_size; y++)
  {
    for (int x = 0; x < lattice_size; x++)
    {
      energy += -spins(y, x)*(spins(y, boundary(x-1)) + spins(boundary(y-1), x));
      magnetization += spins(y, x);
    }
  }
  transition_probabilites.fill(1); // All other (possible) transitions are accepted!
  transition_probabilites(16) = boltzmann(8.0); // simpler than a loop init :)
  transition_probabilites(12) = boltzmann(4.0);
  update_averages(0);
  //Set susceptibility, heatcap
}

int IsingModel::boundary(int r)
{
  // refactor this part!
  int position = r;
  if (r == lattice_size){position = 0;}
  else if(r == -1){position = lattice_size-1;}
  return position;
}

double IsingModel::energy_difference(int x, int y)
{
  return 2*spins(y,x)*(spins(boundary(y-1), x) + spins(boundary(y+1),x)
         + spins(y, boundary(x-1)) + spins(y, boundary(x+1)));
}

double IsingModel::boltzmann(double energy)
{
  return exp(-beta*energy);
}

void IsingModel::update_averages(int cycle)
{
  // update all expectation values
  mean_mag(cycle) =  magnetization;
  mean_abs_mag(cycle) = fabs(magnetization);
  mean_energy(cycle) = energy;
  accepted_flips(cycle) = flips;
}

void IsingModel::finalize()
{
  arma::vec n_values = arma::regspace(1, mc_cycles + 1)*gridsize;
  arma::vec mean_squared_energy;
  arma::vec mean_squared_mag;

  mean_squared_energy = arma::square(mean_energy);
  mean_squared_mag = arma::square(mean_mag);

  mean_mag = arma::cumsum(mean_mag)*1.0/n_values;
  mean_abs_mag = arma::cumsum(mean_abs_mag)*1.0/n_values;
  mean_squared_mag = arma::cumsum(mean_squared_mag)*1.0/n_values;

  all_energies = mean_energy; // save all energies for analysing energy_dist
  mean_energy = arma::cumsum(mean_energy)*1.0/n_values;
  mean_squared_energy = arma::cumsum(mean_squared_energy)*1.0/n_values;

  susceptibility = beta*(arma::square(mean_mag)- mean_squared_mag);
  double factor = (double) beta/temperature;
  heat_capacity = factor*(arma::square(mean_energy)-mean_squared_energy);
}

void IsingModel::metropolis(int seed)
{
  std::random_device rd; // Random seed generator
  std::mt19937 generator(rd() + seed); // Mersenne twister RNG
  std::uniform_int_distribution<> int_dist(0, lattice_size-1); // Uniform integer distribution
  std::uniform_real_distribution<> float_dist(0, 1); // Uniform float distribution

  int x = 0;
  int y = 0;
  float w = 0;
  double current_energy = 0;
  double new_energy = 0;
  double de = 0;
  // Advance simulation one MC cycle
  for(int k = 1; k <= mc_cycles; k++)
  {
    for(int i=0; i < gridsize; i++)
    {
      x = int_dist(generator);
      y = int_dist(generator);

      de = energy_difference(x, y);
      w = transition_probabilites(de + 8); // precalculated transition probability
      if(float_dist(generator) <= w)
      {
        spins(y, x) = -1*spins(y, x);
        energy += (double) de;
        magnetization += 2.0*spins(y, x);
        flips += 1;
      }
    }
    update_averages(k);
  }
  finalize(); // take cumulative sum to find total average
}

void IsingModel::set_ordered_spin_config(std::string config)
{
  if (config == "all_up")
  {
    spins.fill(1);
  }
  else if(config == "all_down")
  {
    spins.fill(-1);
  }
  else if (config == "custom")
  {
    std::cout << "Warning: Custom spin configuration requested!" << std::endl;
  }
  else
  {
    std::cout <<"Error: Spins must be all_up, all_down, or custom."<< std::endl;
    exit(1);
  }
  initialize_parameters(); // reset all quantities.
}

void IsingModel::set_random_spin_config(int seed)
{
    double random_number = 0;
    std::srand(std::time(NULL)+seed);
    for(int y=0; y < lattice_size; y++)
    {
     for(int x=0; x < lattice_size; x++)
     {
       random_number = (double) std::rand()/RAND_MAX;
       spins(y,x) = (random_number < 0.5) ? -1 : 1; // if less than 0.5, down, else up
     }
    }
    initialize_parameters(); // reset all quantities.
}

void IsingModel::set_temperature(double t)
{
  temperature = t;
  initialize_parameters();
}
