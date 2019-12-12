#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"

QuantumDots::QuantumDots(int n, double alpha, double freq, int proc):
   averages(n+1, 5),
   positions(n+1, 6),
   trial_position(6)
    {
      varparam = alpha;
      seed = proc;
      mc_cycles = n;
      omega = freq;
      reset();
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
  double stepsize = (double) sqrt(log(2.0)/(2.0*varparam*omega));

  reset(); // ensure all matrices are emptied

  // initialize positions on unit sphere
  for(int k = 0; k < 6; k++)
  {
    positions(0, k) = 2*float_dist(generator) - 1;
  }

  r = arma::norm(positions.row(0));
  r12 = particle_separation(positions.row(0));

  update_energies(r, r12, 0);
  averages(0, 4) = r12;

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
      update_energies(r, r12, i);
      averages(i, 4) = r12;
    }
    else
    {
      averages.row(i) = averages.row(i-1);
      positions.row(i) = positions.row(i-1);
    }
  }
  finalize();
  //std::cout << (double) accepted_moves/mc_cycles << std::endl;
}

void QuantumDots::reset()
{
  averages.fill(0);
  positions.fill(0);
  trial_position.fill(0);
  accepted_moves = 0;
}

void QuantumDots::finalize()
{
  arma::vec n_values = arma::regspace(1, mc_cycles + 1);

  averages.col(0) = averages.col(2) + averages.col(3); // Etr = Ek + Ep
  averages.col(1) = arma::square(averages.col(0));
  for(int i = 0; i < 5; i++)
  {
    averages.col(i) = arma::cumsum(averages.col(i))/n_values;
  }
}

double QuantumDots::particle_separation(arma::rowvec pos)
{
  return arma::norm(pos.cols(0,2) - pos.cols(3,5));
}

// Psi1 only needs the QuantumDots constructor
Psi1::Psi1(int n, double alpha, double freq, int seed):
          QuantumDots(n, alpha, freq, seed){}

double Psi1::transition_prob(double r, double rp)
{
  return exp(-varparam*omega*(rp*rp - r*r));
}

void Psi1::update_energies(double r, double r12, int cycle)
{
  double factor = 0.5*omega*omega*r*r;
  double ek = -factor*varparam*varparam + 3.0*varparam*omega;
  double ep = factor + 1.0/r12;
  averages(cycle, 2) = ek;
  averages(cycle, 3) = ep;
}

Psi2::Psi2(int n, double alpha, double beta , double freq, int seed):
      QuantumDots(n, alpha, freq, seed)
{
  a = alpha;
  varparam = beta;
}

double Psi2::transition_prob(double r, double rp)
{
  double r12p = particle_separation(trial_position);
  double frac = (double) r12p/(1+varparam*r12p) - (double) r12/(1+varparam*r12);
  return exp(-a*omega*(rp*rp - r*r) + frac);
}

void Psi2::update_energies(double r, double r12, int cycle)
{
  double temp = 0.5*omega*omega*r*r;
  double ek1 = -temp*a*a + 3.0*a*omega;
  double factor = 1.0 + varparam*r12;
  double frac = 1.0/(2.0*factor*factor);
  double ek = ek1 + frac*(a*omega*r12 - frac - 2.0/r12 + varparam*2.0/factor);
  double ep = temp + 1.0/r12;
  averages(cycle, 2) = ek;
  averages(cycle, 3) = ep;
}

/*
double Psi2::local_energy(double r, double r12)
{
  double factor = 1.0+varparam*r12;
  double frac = 1.0/(2.0*factor*factor);
  double e1 = 0.5*omega*omega*r*r*(1-a*a) + 3*a*omega + 1.0/r12;
  return e1 + frac*(a*omega*r12 - frac - 2.0/r12 + varparam*2.0/factor);
}



double Psi1::local_ep(double r, double r12)
{
  return 0.5*omega*omega*r*r*(1-varparam*varparam) + 3*varparam*omega + 1.0/r12;
}
*/
