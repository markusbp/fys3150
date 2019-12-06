#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"

QuantumDots::QuantumDots(int n, double alpha, double freq, int proc):
   averages(n+1, 4),
   positions(n+1, 6)
    {
      varparam = alpha;
      seed = proc;
      mc_cycles = n;
      omega = freq;
      stepsize = (double) sqrt(log(2.0)/(2.0*varparam*omega));
      averages.fill(0);
      positions.fill(0);
    }

void QuantumDots::metropolis()
{
  std::random_device rd; // Random seed generator, extra seed can be added
  std::mt19937 generator((double) rd()/(1.0 + seed)); // Mersenne twister RNG
  std::uniform_real_distribution<> float_dist(0, 1); // Uniform float distribution

  double w = 0;
  double s = 0;
  double r = 0;
  double rp = 0;
  double r12 = 0;
  double accepted_moves = 0;
  arma::rowvec trial_position(6, arma::fill::zeros); //x1, x2, y1, y2, z1, z2
  // initialize positions
  for(int k = 0; k < 6; k++)
  {
    positions(0, k) = 2*float_dist(generator) - 1;
  }

  r = arma::norm(positions.row(0));
  r12 = particle_separation(positions.row(0));

  averages(0, 0) = local_energy(r, r12);
  averages(0, 1) = r12;

  for(int i = 1; i <= mc_cycles; i++)
  {
    for(int j = 0; j < 6; j++)
    {
      trial_position(j) = positions(i-1, j)+(2*float_dist(generator)- 1)*stepsize;
    }

    rp = arma::norm(trial_position);
    w = transition_prob(r, rp);
    s = float_dist(generator);

    if (s < w)
    {
      accepted_moves += 1.0;
      positions.row(i) = trial_position;
      r = rp;
      r12 = particle_separation(trial_position);
      averages(i, 0) = local_energy(r, r12);
      averages(i, 1) = r12;
    }
    else
    {
      averages(i) = averages(i-1);
      positions.row(i) = positions.row(i-1);
    }
  }
  finalize();
  std::cout << (double) accepted_moves/mc_cycles << std::endl;
}

void QuantumDots::finalize()
{
  arma::vec n_values = arma::regspace(1, mc_cycles + 1);
  for(int i = 0; i < 4; i++)
  {
    averages.col(i) = arma::cumsum(averages.col(i))/n_values;
  }
}

double QuantumDots::particle_separation(arma::rowvec pos)
{
  return arma::norm(pos.cols(0,2) - pos.cols(3,5));
}

double QuantumDots::transition_prob(double r, double rp)
{
  return exp(-varparam*omega*(rp*rp - r*r));
}

double QuantumDots::local_energy(double r, double r12)
{
  return 0.5*omega*omega*r*r*(1-varparam*varparam) + 3*varparam*omega + 1.0/r12;
}

void QuantumDots::grid_search(int steps, double step)
{
  double min_en  = 1e6;
  double min_par = 0;
  averages.fill(0);
  positions.fill(0);
  for(int i = 0; i<steps; i++)
  {
    stepsize = (double) sqrt(log(2.0)/(2.0*varparam*omega));
    metropolis();

    if(averages(mc_cycles, 0) < min_en)
    {
      min_en = averages(mc_cycles,0);
      min_par = varparam;
    }

    averages.fill(0);
    positions.fill(0);
    // reset
    varparam += step;
  }
  std::cout << min_par << " En: "<< min_en << std::endl;

}
