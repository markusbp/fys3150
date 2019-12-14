#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>
#include "quantum_dots.hpp"

QuantumDots::QuantumDots(int n, double alpha, double freq, int proc):
   averages(n+1, 7),
   positions(n+1, 6),
   trial_position(6)
    { // constructor for abstract class
      varparam = alpha; // desired variational parameter
      seed = proc; // most often used for process number with MPI
      mc_cycles = n;
      omega = freq; // HO frequency
      reset(); // fill all matrices with zero
    }

void QuantumDots::metropolis()
{ // run MC simulation based on metropolis algorithm
  std::random_device rd; // Random seed generator, extra seed can be added
  std::mt19937 generator((double) rd()/(1.0 + seed)); // Mersenne twister RNG
  std::uniform_real_distribution<> float_dist(0, 1); // Uniform float distribution

  double w = 0; // transition probability
  double s = 0; // random number used to reject or accept moves
  double r = 0; // length of total particle position(s) vector(s)
  double rp = 0; // length of proposed position vector
  double stepsize = find_stepsize(); // maximum stepsize

  reset(); // ensure all matrices are emptied

  // initialize positions on unit sphere
  for(int k = 0; k < 6; k++)
  {
    positions(0, k) = 2*float_dist(generator) - 1; // random step between -1 & 1
  }

  r = arma::norm(positions.row(0)); // initial position
  r12 = particle_separation(positions.row(0)); // initial separation

  update_energies(r, r12, 0); // initial energies
  averages(0, 4) = r12; // save separation

  for(int i = 1; i <= mc_cycles; i++)
  {
    for(int j = 0; j < 6; j++)
    { // propose new random position
      trial_position(j) = positions(i-1, j)+(2*float_dist(generator)- 1)*stepsize;
    }

    rp = arma::norm(trial_position);
    w = transition_prob(r, rp); // find probability of move
    s = float_dist(generator);

    if (s < w)
    { // move is accepted: update all quantities
      accepted_moves += 1.0;
      positions.row(i) = trial_position;
      r = rp;
      r12 = particle_separation(trial_position);
      update_energies(r, r12, i);
      averages(i, 4) = r12;
    }
    else
    { // move is rejected, stay at previous parameters
      averages.row(i) = averages.row(i-1);
      positions.row(i) = positions.row(i-1);
    }
  }
  finalize(); // calculate averages
}

void QuantumDots::reset()
{ // clear all
  averages.fill(0);
  positions.fill(0);
  trial_position.fill(0);
  accepted_moves = 0;
}

void QuantumDots::finalize()
{ // calculate all averages
  arma::vec n_values = arma::regspace(1, mc_cycles + 1);
  arma::vec separation = 1.0/averages.col(4);
  // Non-interacting energy: Ek + Ep_noninteracting
  averages.col(5) = averages.col(2) + averages.col(3); // non-interacting
  averages.col(6) = averages.col(3); // save non-interacting potential
  // full energy = Ek + Ep_noninteracting + 1.0/r12
  averages.col(0) = averages.col(5) + separation; // full energy
  averages.col(3) = averages.col(3) + separation; // ep with interaction
  averages.col(1) = arma::square(averages.col(0)); // squared energy
  for(int i = 0; i < 7; i++)
  { // cumsum/n_i gives averages
    averages.col(i) = arma::cumsum(averages.col(i))/n_values;
  }
}

double QuantumDots::particle_separation(arma::rowvec pos)
{ // find particle separation sqrt((x1-x2)^2+(y1-y2)^2+...
  return arma::norm(pos.cols(0,2) - pos.cols(3,5));
}

// Psi1 only needs the QuantumDots constructor
Psi1::Psi1(int n, double alpha, double freq, int seed):
          QuantumDots(n, alpha, freq, seed){}

double Psi1::transition_prob(double r, double rp)
{ // Psi1 transition probability
  return exp(-varparam*omega*(rp*rp - r*r));
}

void Psi1::update_energies(double r, double r12, int cycle)
{ // psi1 energies
  double ep_nonint = 0.5*omega*omega*r*r; // noninteracting harmonic oscillator
  double ek = -ep_nonint*varparam*varparam + 3.0*varparam*omega;
  averages(cycle, 2) = ek; // kinetic energy
  averages(cycle, 3) = ep_nonint; // ep, noninteracting, interaction added later
}

double Psi1::find_stepsize()
{ // psi1 step size
  return (double) sqrt(log(2.0)/(2.0*varparam*omega));
}

Psi2::Psi2(int n, double alpha, double beta, double freq, int seed):
      QuantumDots(n, alpha, freq, seed)
{ // Psi2 constructor
  a = alpha;
  varparam = beta;
}

double Psi2::transition_prob(double r, double rp)
{ // Psi2 transition probability
  double r12p = particle_separation(trial_position);
  double frac = (double) r12p/(1+varparam*r12p) - (double) r12/(1+varparam*r12);
  return exp(-a*omega*(rp*rp - r*r) + frac);
}

void Psi2::update_energies(double r, double r12, int cycle)
{ // update potential and kinetic energy, psi2
  double ep_nonint = 0.5*omega*omega*r*r; // noninteracting potential energy
  double ek1 = -ep_nonint*a*a + 3.0*a*omega;
  double factor = 1.0 + varparam*r12;
  double frac = 1.0/(2.0*factor*factor);
  double ek = ek1 + frac*(a*omega*r12 - frac - 2.0/r12 + varparam*2.0/factor);
  averages(cycle, 2) = ek; // kinetic energy
  averages(cycle, 3) = ep_nonint; // ep, noninteracting, interaction added later
}

double Psi2::find_stepsize()
{ // psi2 stepsize
  return (double) sqrt(log(2.0)/(2.0*a*omega));
}

void Psi2::set_omega(double freq)
{ // psi2 setter method for frequency
  omega = freq;
}
