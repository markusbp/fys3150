#ifndef ISING_HPP
#define ISING_HPP

class IsingModel
{
  private:
    double kb; // Boltzmann constant
    double beta; // inverse temperature
    double boltzmann(double energy); // boltzmann factor
    double energy_difference(int x, int y); // calculate energy difference
    int boundary(int r); // periodic boundary conditions
    void update_averages(int cycle);
    void initialize_parameters(); // energy and magnetization
    void finalize(); // calculate expected value

  public:
    int mc_cycles;  // Monte Carlo cycles
    int lattice_size;  // lxl grid
    int gridsize;
    double temperature; // temperature
    double energy;
    double magnetization; // mean magnetization
    double flips;
    arma::imat spins; // init spin matrix
    arma::vec transition_probabilites; // probabilities, save calculations
    arma::vec susceptibility;
    arma::vec heat_capacity;
    arma::vec mean_mag;
    arma::vec mean_abs_mag;
    arma::vec mean_energy;
    arma::vec all_energies;
    arma::vec accepted_flips;
    void metropolis(int seed); // method for running step of Metropolis algorithm
    void set_ordered_spin_config(std::string config); // helper methods for setting
    void set_random_spin_config(int seed); // spin config
    void set_temperature(double t); // or temperature
    IsingModel(int l, int n, double t); // constructor
};

#endif
