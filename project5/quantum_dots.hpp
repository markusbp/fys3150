#ifndef QDOTS_HPP
#define QDOTS_HPP

class QuantumDots
{
  protected:
    int seed;
    double omega;
    void finalize();
    double particle_separation(arma::rowvec pos);
    virtual double transition_prob(double r1, double r2);
    virtual double local_energy(double r1, double r12);
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
    double local_energy(double r1, double r12);
  public:
    Psi1(int n, double alpha, double freq, int seed);
};

class Psi2: public QuantumDots
{
  protected:
    double a;
    double transition_prob(double r1, double r2);
    double local_energy(double r1, double r12);
  public:
    Psi2(int n, double alpha, double beta, double freq, int seed);
};

#endif
