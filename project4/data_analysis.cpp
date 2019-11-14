#include <armadillo>
#include "ising_model.hpp"


// file for creating all data in task.
void save_all(IsingModel model, std::string filename, std::string t)
{
  model.mean_abs_mag.save(filename + "_mean_abs_mag" + "_t=" + t, arma::csv_ascii);
  model.mean_mag.save(filename + "_mean_mag"+ "_t=" + t, arma::csv_ascii);
  model.mean_energy.save(filename + "_mean_energy"+ "_t=" + t, arma::csv_ascii);
  model.heat_capacity.save(filename + "_heatcap"+ "_t=" + t, arma::csv_ascii);
  model.susceptibility.save(filename + "_susceptibility"+ "_t=" + t, arma::csv_ascii);
  model.accepted_flips.save(filename + "_accepted_flips"+ "_t=" + t, arma::csv_ascii);
}

void save_en_mag(IsingModel model, std::string filename, std::string t)
{
  model.mean_abs_mag.save(filename + "_mean_abs_mag" + "_t=" + t, arma::csv_ascii);
  model.mean_energy.save(filename + "_mean_energy"+ "_t=" + t, arma::csv_ascii);
}

int main(int argc, char *argv[])
{
  // 2x2, task 4b)

  int l = 2;
  int n = 1e6;
  double t = 1.0;
  /*
  IsingModel model(l, n, t); //  l, n, temp
  model.metropolis(0);
  save_all(model, "./results/2x2/ordered", "1.0");

  model.set_random_spin_config(0);
  model.metropolis(0);
  save_all(model, "./results/2x2/random", "1.0");
  */
  // 20x20
  l = 20;
  n = 1e6;
  t = 2.3;
  IsingModel model20(l, n, t); //  l, n, temp

  model20.metropolis(0);
  model20.all_energies.save("./results/20x20/all_energies", arma::csv_ascii);

  /*
  save_all(model20, "./results/20x20/ordered", "1.0");

  model20.set_random_spin_config(0);
  model20.metropolis(0);
  save_all(model20, "./results/20x20/random", "1.0");

  model20.set_temperature(2.4); // set temp and reset values

  model20.metropolis(0);
  save_all(model20, "./results/20x20/ordered", "2.4");

  model20.set_random_spin_config(0);
  model20.metropolis(0);
  save_all(model20, "./results/20x20/random", "2.4");
  */
  return 0;
}
