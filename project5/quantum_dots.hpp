#ifndef QDOTS_HPP
#define QDOTS_HPP

class QuantumDots
{
  private:
    int mc_cycles;
    int seed;
    double omega;
    void initialize();
    void finalize();
    double particle_separation(arma::rowvec pos);
  public:
    double a;
    double b;
    double stepsize;
    arma::mat averages;
    arma::mat positions;
    void metropolis();
    void grid_search(int steps, double step);
    double transition_prob(double r1, double r2);
    double local_energy(double r1, double r12);
    QuantumDots(int n, double alpha, double beta, double freq, int seed);
};

#endif
