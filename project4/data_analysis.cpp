#include <armadillo>
#include "ising_model.hpp"

// File for generating data for tasks b-c-d
// Uses Ising Model class
// All results printed to file, no real printout to show

// usage:
// ./analysis task -temperature
// -temperature is only required when task = 'correlation'
// task is either 'b', for task b, 'cd' for tasks c and d, or 'correlation'
// correlation simulates and saves spin config for given optional -temperature

// file for creating all data in task.
void save_params(IsingModel model, std::string filename, std::string t)
{ // helper function for saving important parameters
  t = t.substr(0, 4);
  model.mean_abs_mag.save(filename + "_mean_abs_mag" + "_t=" + t, arma::csv_ascii);
  model.mean_mag.save(filename + "_mean_mag"+ "_t=" + t, arma::csv_ascii);
  model.mean_energy.save(filename + "_mean_energy"+ "_t=" + t, arma::csv_ascii);
  model.heat_capacity.save(filename + "_heatcap"+ "_t=" + t, arma::csv_ascii);
  model.susceptibility.save(filename + "_susceptibility"+ "_t=" + t, arma::csv_ascii);
}

void save_en_mag(IsingModel model, std::string filename, std::string t)
{ // helper function for saving important parameters
  t = t.substr(0, 4);
  model.mean_abs_mag.save(filename + "_mean_abs_mag" + "_t=" + t, arma::csv_ascii);
  model.mean_energy.save(filename + "_mean_energy"+ "_t=" + t, arma::csv_ascii);
  model.accepted_flips.save(filename + "_accepted_flips" + "_t=" + t, arma::csv_ascii);
  model.all_energies.save(filename + "_all_energies" + "_t=" + t, arma::csv_ascii);
}

int main(int argc, char *argv[])
{
  std::string task = argv[1];
  std::string temp;

  int l;
  int n;
  double t;

  if(task == "b")
  {  // 2x2 lattice, task 4b)
    l = 2; // grid size, 2x2
    n = 1e5; // number of monte carlo cycles
    t = 1.0; // temperature, first run
    temp = std::to_string(t);
    IsingModel model(l, n, t); //  l, n, temp
    model.metropolis(0); // run simulation with only one process.
    save_params(model, "./results/2x2/", temp);
  }
  if(task == "cd")
  {
    // Runs Ising Model and creates data for tasks 4c and d,
    l = 20; // 20x20 lattice
    n = 1e6; // mc cycles
    t = 1.0; // initial temp
    temp = std::to_string(t);

    IsingModel model(l, n, t); // instantiate model
    model.metropolis(0); // run metropolis algo for n mc cycles
    save_en_mag(model, "./results/20x20/ordered", temp); // save relevant params

    model.set_random_spin_config(0); // Do the same for random spin config; reset
    model.metropolis(0);
    save_en_mag(model, "./results/20x20/random", temp);

    t = 2.4;
    model.set_temperature(t); // Perform same calculations for temperature 2.4
    temp = std::to_string(t);

    model.set_ordered_spin_config("all_up");
    model.metropolis(0);
    save_en_mag(model, "./results/20x20/ordered", temp);

    model.set_random_spin_config(0);
    model.metropolis(0);
    save_en_mag(model, "./results/20x20/random", temp);
  }
  if (task == "correlation")
  {
    l = 100;
    n = 1e5;
    t = atof(argv[2]);
    IsingModel model(l, n, t); // instantiate model
    model.set_temperature(t);
    model.set_random_spin_config(0);
    model.metropolis(0); //simulate and save only spin configuration
    model.spins.save("./results/spin/spins_t=" + std::to_string(t), arma::csv_ascii);
  }
  return 0;
}
