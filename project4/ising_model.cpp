#include <iostream>
#include <random>
#include <cmath>
#include <armadillo>

class IsingModel
{
  private:
    double kb; // Boltzmann constant, m2kg2s-2K-1
    double beta;
    double random_lattice_site(){return 0;} // fetch random site on lattice
    int calculate_energy(int x, int y);
    double boltzmann(double energy);
    void update_expectation_values(int x, int y);

  public:
    int n;  // Monte Carlo cycles
    int l;  // lxl grid
    double t;
    double energy;
    double magnetization; // mean (absolute) magnetization
    arma::imat spins; // init spin matrix

    void metropolis();

    IsingModel(int mc_cycles, int lattice_size, double temperature):
              spins(lattice_size + 2, lattice_size + 2) // Constructor
    {
      kb = 1.3806e-23;
      t = temperature;
      beta = (double) 1.0/(kb*t);
      n = mc_cycles;
      l = lattice_size;
      spins.fill(1); // Default configuration
    }
};

void IsingModel::metropolis()
{
  std::random_device rd; // Random seed generator
  std::mt19937 generator(rd()); // Mersenne twister RNG
  std::uniform_int_distribution<> int_dist(1, l-1); // Uniform integer distribution
  std::uniform_real_distribution<> float_dist(0, 1); // Uniform integer distribution

  int x = 0;
  int y = 0;
  float w = 0;
  double new_energy = 0;
  double de = 0;

  for(int i=0; i < l*l*n; i++)
  {
    x = int_dist(generator);
    y = int_dist(generator);
    new_energy = calculate_energy(x, y);

    de = new_energy - energy;

    if(de < 0)
    {
      spins(x,y) = -spins(x, y);
    }
    else if(de > 0)
    {
      w = boltzmann(de);
      if(float_dist(generator) <= w)
      {
        spins(x, y) = -spins(x, y);
      }
    }
    update_expectation_values(x, y);
  }
  // REMEMBER TO CHECK BOUNDARIES!!
  // flip random spin, if energetically favourable, flip, else only possibly flip, according to Boltzmann distribution
}

int IsingModel::calculate_energy(int x, int y)
{
  return -2*spins(x,y)*(spins(x-1, y) + spins(x+1,y) + spins(x, y-1) + spins(x, y+1));
}

double IsingModel::boltzmann(double delta_e)
{
  return exp(-beta*delta_e);
}

void IsingModel::update_expectation_values(int x, int y)
{
  // update all expectation values
  magnetization += 2*spins(x, y);

}

int main()
{
  IsingModel model(10, 40, 2);
  //std::cout << model.spins << std::endl;
  model.metropolis();
  std::cout << model.spins << std::endl;
  return 0;
}
