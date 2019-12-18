#ifndef QDOTS_HPP
#define QDOTS_HPP

// Class implementation of 2 electron quantum dot system, in Harmonic Oscillator
// potential, with different trial wavefunctions, for use with variational
// Monte Carlo simulation.

class QuantumDots
{
  // Abstract class, made to be inherited for wave function-specific class
  protected:
    int seed; // seed for each process
    double r12; // particle separation
    double omega; // Harmonic oscillator frequency
    void finalize(); // compute averages
    double particle_separation(arma::rowvec pos); // compute r12
    virtual double transition_prob(double r1, double r2) = 0; // transition probability
    virtual double find_stepsize() = 0; // determine stepsize for a given model
    virtual void update_energies(double r1, double r12, int cycle) = 0;
    arma::rowvec trial_position; //x1, x2, y1, y2, z1, z2 (coordinates)
  public:
    int mc_cycles; // number of simulation cycles
    double varparam; // variational parameter
    double accepted_moves; // number of accepted mc moves
    arma::mat averages; // all interesting quantities
    arma::mat positions; // all particle positions
    void reset(); // clear all
    void metropolis(); // metropolis algorithm MC simulation
    QuantumDots(int n, double alpha, double freq, int seed); // constructor
};

class Psi1: public QuantumDots
{ // same methods as for QuantumDots, because inheritance :)
  protected:
    double transition_prob(double r1, double r2);
    double find_stepsize();
    void update_energies(double r1, double r12, int cycle);
  public:
    Psi1(int n, double alpha, double freq, int seed);
    void set_omega(double freq); // change omega without redoing all classes :)
};

class Psi2: public QuantumDots
{
  protected:
    double transition_prob(double r1, double r2);
    double find_stepsize();
    void update_energies(double r1, double r12, int cycle);
  public:
    double a; // a is public do to better grid search
    void set_omega(double freq); // change omega without redoing all classes :)
    Psi2(int n, double alpha, double beta, double freq, int seed);
};

#endif
