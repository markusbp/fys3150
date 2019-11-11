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
    double boltzmann(double energy); // boltzmann factor
    int boundary(int r);
    int calculate_energy(int x, int y, int sign); // calculate energy difference
    void update_averages(int x, int y);
    void find_initial_parameters(); // energy and magnetization

  public:
    int n;  // Monte Carlo cycles
    int l;  // lxl grid
    double t;
    double energy;
    double magnetization; // mean (absolute) magnetization
    double abs_magnetization;
    arma::imat spins; // init spin matrix
    arma::vec transition_probabilites; // probabilities, save calculations
    void metropolis();

    IsingModel(int mc_cycles, int lattice_size, double temperature):
              spins(lattice_size, lattice_size), transition_probabilites(17) // Constructor
    {
      kb =  1;//1.3806e-23;
      t = temperature;
      beta = (double) 1.0/(kb*t);
      n = mc_cycles;
      l = lattice_size;
      spins.fill(1); // Default configuration
      energy = 0;
      magnetization = 0;
      abs_magnetization = 0;
      find_initial_parameters();
    }
};

void IsingModel::find_initial_parameters()
{
  for(int y = 0; y < l; y++)
  {
    for (int x = 0; x < l; x++)
    {
      energy += -spins(y, x)*(spins(y, boundary(x-1)) + spins(boundary(y-1), x));
      magnetization += spins(y, x);
    }
  }
  transition_probabilites.fill(1); // All other transitions are accepted!
  transition_probabilites(16) = boltzmann(8.0); // simpler than a loop init :)
  transition_probabilites(12) = boltzmann(4.0);
}

int IsingModel::boundary(int r)
{
  // refactor this part!
  int position = r;
  if (r == l){position = 0;}
  else if(r == -1){position = l-1;}
  return position;
}

int IsingModel::calculate_energy(int x, int y, int sign)
{
  return -sign*spins(y,x)*(spins(boundary(y-1), x) + spins(boundary(y+1),x)
         + spins(y, boundary(x-1)) + spins(y, boundary(x+1)));
}

double IsingModel::boltzmann(double delta_e)
{
  return exp(-beta*delta_e);
}

void IsingModel::update_averages(int x, int y)
{
  // update all expectation values
  magnetization += 2*spins(y, x);
  abs_magnetization += fabs(magnetization);
  std::cout << (double)magnetization/(n*l*l) << std::endl;;
}

void IsingModel::metropolis()
{
  std::random_device rd; // Random seed generator
  std::mt19937 generator(rd()); // Mersenne twister RNG
  std::uniform_int_distribution<> int_dist(0, l-1); // Uniform integer distribution
  std::uniform_real_distribution<> float_dist(0, 1); // Uniform integer distribution

  int x = 0;
  int y = 0;
  float w = 0;
  double current_energy = 0;
  double new_energy = 0;
  double de = 0;

  for(int i=0; i < l*l*n; i++)
  {
    x = int_dist(generator);
    y = int_dist(generator);

    current_energy = calculate_energy(y, x, 1);
    // do a flip
    new_energy = calculate_energy(y, x, -1); // calculate new energy
    de = new_energy - current_energy; // find energy change
    w = transition_probabilites(de + 8);
    if(float_dist(generator) < w)
    {
      spins(y, x) = -spins(y, x);
      energy += de;
      update_averages(y, x);
    }
  }
  // flip random spin, if energetically favourable, flip, else only possibly flip, according to Boltzmann distribution
}

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]);
  int l = atoi(argv[2]);
  float t = atof(argv[3]);
  IsingModel model(n, l, t); // n, l, temp
  model.metropolis();

  return 0;
}
