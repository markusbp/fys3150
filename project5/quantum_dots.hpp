#ifndef QDOTS_HPP
#define QDOTS_HPP

class QuantumDots
{
  protected:
    int seed;
    double r12;
    double omega;
    double accepted_moves;
    void finalize();
    double particle_separation(arma::rowvec pos);
    virtual double transition_prob(double r1, double r2) = 0;
    virtual void update_energies(double r1, double r12, int cycle) = 0;
    arma::rowvec trial_position; //x1, x2, y1, y2, z1, z2
  public:
    int mc_cycles;
    double varparam;
    arma::mat averages;
    arma::mat positions;
    void reset();
    void metropolis();
    QuantumDots(int n, double alpha, double freq, int seed);
};

class Psi1: public QuantumDots
{
  protected:
    double transition_prob(double r1, double r2);
    void update_energies(double r1, double r12, int cycle);
  public:
    Psi1(int n, double alpha, double freq, int seed);
};

class Psi2: public QuantumDots
{
  protected:
    double transition_prob(double r1, double r2);
    void update_energies(double r1, double r12, int cycle);
  public:
    double a; // a is public do to better grid search
    void set_omega(double freq); // change omega without redoing all classes :)
    Psi2(int n, double alpha, double beta, double freq, int seed);
};

#endif
